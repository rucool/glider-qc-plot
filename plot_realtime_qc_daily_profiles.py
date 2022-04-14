#!/usr/bin/env python

"""
Author: Lori Garzio on 2/4/2022
Last modified: 2/24/2022
Plot realtime glider science profiles from yesterday, values flagged by QC variables are highlighted. Plots CTD and
dissolved oxygen variables (if available).
"""

import argparse
import sys
import numpy as np
import pandas as pd
import datetime as dt
import os
import matplotlib.pyplot as plt
import functions.common as cf
plt.rcParams.update({'font.size': 12})


def define_markers(qc_varname):
    markers = dict(climatology=dict(m='v', s=60, alpha=1),
                   hysteresis=dict(m='s', s=40, alpha=1),
                   flat_line=dict(m='^', s=60, alpha=1),
                   gross_range=dict(m='D', s=40, alpha=1),
                   rate_of_change=dict(m='X', s=80, alpha=1),
                   spike=dict(m='*', s=100, alpha=1),
                   pressure=dict(m='s', s=40, alpha=1),
                   summary=dict(m='o', s=100, alpha=.2)
                   )
    mkey = [key for key in markers.keys() if key in qc_varname][0]
    return markers[mkey]


def main(args):
    sdir = args.save_dir
    test = args.test
    profile_groups = args.profile_groups

    if test:
        erddap_server = 'http://slocum-test.marine.rutgers.edu//erddap'
    else:
        erddap_server = 'http://slocum-data.marine.rutgers.edu//erddap'

    for deployment in args.deployments:
        year = dt.datetime.strptime(deployment.split('-')[-1], '%Y%m%dT%H%M').year
        glider_id = '{}-profile-sci-rt'.format(deployment)

        # get yesterday's date
        today = dt.date.today()
        # today = dt.date.today() - dt.timedelta(days=340)  # for debugging
        t0 = today - dt.timedelta(days=1)
        t0_savestr = t0.strftime('%Y%m%d')
        t0str = t0.strftime('%Y-%m-%dT00:00:00Z')
        t1str = today.strftime('%Y-%m-%dT00:00:00Z')

        constraints = dict({'time>=': t0str, 'time<=': t1str})

        kwargs = dict()
        kwargs['constraints'] = constraints
        ds = cf.get_erddap_dataset(erddap_server, glider_id, **kwargs)
        ds = ds.swap_dims({'obs': 'time'})
        ds = ds.sortby(ds.time)

        savedir_dt = os.path.join(sdir, str(year), glider_id, t0_savestr)
        os.makedirs(savedir_dt, exist_ok=True)

        profiletimes = np.unique(ds.profile_time.values)
        num_profiles = len(profiletimes)

        plot_sections = np.arange(0, num_profiles, profile_groups)
        plot_sections = np.append(plot_sections, num_profiles)

        varlist = []
        for il in ['ctd', 'do']:
            varlist.append(cf.define_instrument_variables(il))
        varlist = cf.flatten(varlist)

        var_list = set(ds.data_vars).intersection(varlist)

        flag_defs = cf.flag_defs()

        for ps_idx, ps in enumerate(plot_sections):
            if np.logical_and(ps == num_profiles, len(plot_sections) > 2):
                continue
            if ps_idx > 0:
                if ps_idx == 1:
                    ii = 0
                else:
                    ii = plot_sections[ps_idx - 1] + 1
                ptimes = profiletimes[ii:ps]
                ptimes_idx = np.where(np.logical_and(ds.profile_time >= ptimes[0], ds.profile_time <= ptimes[-1]))[0]
                time0 = np.nanmin(ds.time.values[ptimes_idx])
                time1 = np.nanmax(ds.time.values[ptimes_idx])
                dss = ds.sel(time=slice(time0, time1))
                t0str = pd.to_datetime(np.nanmin(dss.profile_time.values)).strftime('%Y-%m-%dT%H:%M')
                t1str = pd.to_datetime(np.nanmax(dss.profile_time.values)).strftime('%Y-%m-%dT%H:%M')
                t0save = pd.to_datetime(np.nanmin(dss.profile_time.values)).strftime('%Y%m%dT%H%M')
                t1save = pd.to_datetime(np.nanmax(dss.profile_time.values)).strftime('%Y%m%dT%H%M')
                for cv in var_list:
                    save_filename = f'{cv}_qc_{t0save}-{t1save}.png'
                    try:
                        data = dss[cv]
                    except KeyError:
                        continue
                    pressure = dss.pressure
                    pressure_interp = cf.interpolate_pressure(pressure)
                    fig, ax = plt.subplots(figsize=(8, 10))

                    # iterate through each profile and plot the profile lines
                    for pt in ptimes:
                        pt_idx = np.where(dss.profile_time.values == pt)[0]
                        non_nans = np.where(np.invert(np.isnan(data[pt_idx])))[0]
                        ax.plot(data[pt_idx][non_nans], pressure_interp[pt_idx][non_nans], color='gray')  # plot lines

                    # add points
                    ax.scatter(data, pressure_interp, color='gray', s=20, zorder=5)

                    # find the qc variables (don't plot summary)
                    qc_vars = [x for x in ds.data_vars if np.logical_and(f'{cv}_' in x, '_summary_' not in x)]
                    if cv in ['salinity', 'density']:
                        qc_vars.append('conductivity_hysteresis_test')
                        qc_vars.append('temperature_hysteresis_test')

                    for qi, qv in enumerate(qc_vars):
                        try:
                            flag_vals = dss[qv].values
                        except KeyError:
                            continue
                        for fd, info in flag_defs.items():
                            qc_idx = np.where(flag_vals == info['value'])[0]
                            if len(qc_idx) > 0:
                                m_defs = define_markers(qv)
                                ax.scatter(data[qc_idx], pressure_interp[qc_idx], color=info['color'], s=m_defs['s'],
                                           marker=m_defs['m'], edgecolor='k', alpha=m_defs['alpha'],
                                           label=f'{qv}-{fd}', zorder=10)

                    # add legend if necessary
                    handles, labels = plt.gca().get_legend_handles_labels()
                    by_label = dict(zip(labels, handles))
                    if len(handles) > 0:
                        ax.legend(by_label.values(), by_label.keys(), loc='best')

                    ax.invert_yaxis()
                    ax.set_ylabel('Pressure (dbar)')
                    ax.ticklabel_format(useOffset=False)  # don't use scientific notation for ticks
                    ax.set_xlabel(f'{cv} ({data.units})')
                    ttl = f'{deployment} {t0str} to {t1str}'
                    ax.set_title(ttl)

                    sfile = os.path.join(savedir_dt, save_filename)
                    plt.savefig(sfile, dpi=300)
                    plt.close()


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description='Plot real time glider QC data',
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('deployments',
                            nargs='+',
                            help='Glider deployment name(s) formatted as glider-YYYYmmddTHHMM')

    arg_parser.add_argument('-s', '--save_dir',
                            dest='save_dir',
                            type=str,
                            help='Full file path to save directory')

    arg_parser.add_argument('-pg', '--profile_groups',
                            help='number of profiles to group together in 1 plot',
                            type=int,
                            default=10)

    arg_parser.add_argument('-test', '--test',
                            help='Download data from the slocum-test ERDDAP server.',
                            action='store_true')

    parsed_args = arg_parser.parse_args()
    sys.exit(main(parsed_args))
