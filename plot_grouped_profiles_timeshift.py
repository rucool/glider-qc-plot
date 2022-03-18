#!/usr/bin/env python

"""
Author: lgarzio on 3/16/2022
Last modified: lgarzio on 3/16/2022
Plot groups of profiles for time-shifted variables, before and after shift.
Data are plotted from sci-profile datasets downloaded from RUCOOL's glider ERDDAP server using download_dataset.py
"""

import numpy as np
import pandas as pd
import xarray as xr
import os
import matplotlib.pyplot as plt
import functions.common as cf
plt.rcParams.update({'font.size': 12})


def main(ncf, sdir, nprof):
    fname = ncf.split('/')[-1].split('.nc')[0]
    deploy = f'{fname.split("-")[0]}-{fname.split("-")[1]}'

    ds = xr.open_dataset(ncf)
    ds = ds.swap_dims({'obs': 'time'})
    ds = ds.sortby(ds.time)

    savedir = os.path.join(sdir, deploy, fname, f'profiles_timeshift_group{nprof}')
    os.makedirs(savedir, exist_ok=True)

    profiletimes = np.unique(ds.profile_time.values)

    plot_sections = np.arange(0, len(profiletimes), nprof)
    plot_sections = np.append(plot_sections, len(profiletimes))

    shifted_vars = [v for v in list(ds.data_vars) if '_shifted' in v]

    for ps_idx, ps in enumerate(plot_sections):
        if ps_idx > 0:
            if ps_idx == 1:
                ii = 0
            else:
                ii = plot_sections[ps_idx - 1] + 1
            ptimes = profiletimes[ii:ps]
            try:
                ptimes_idx = np.where(np.logical_and(ds.profile_time >= ptimes[0], ds.profile_time <= ptimes[-1]))[0]
            except IndexError:
                continue
            time0 = np.nanmin(ds.time.values[ptimes_idx])
            time1 = np.nanmax(ds.time.values[ptimes_idx])
            dss = ds.sel(time=slice(time0, time1))
            t0str = pd.to_datetime(np.nanmin(dss.profile_time.values)).strftime('%Y-%m-%dT%H:%M')
            t1str = pd.to_datetime(np.nanmax(dss.profile_time.values)).strftime('%Y-%m-%dT%H:%M')
            t0save = pd.to_datetime(np.nanmin(dss.profile_time.values)).strftime('%Y%m%dT%H%M')
            t1save = pd.to_datetime(np.nanmax(dss.profile_time.values)).strftime('%Y%m%dT%H%M')
            for sv in shifted_vars:
                save_filename = f'{sv}_{t0save}-{t1save}.png'

                data_shifted = dss[sv]
                data = dss[sv.split('_shifted')[0]]

                # in some cases, ERDDAP doesn't set the metadata/fill values, so get rid of any possible fill values
                data_shifted[data_shifted > 10000] = np.nan
                data[data > 10000] = np.nan

                # convert 0.0 pH voltage values to nan
                if 'ph_ref_voltage' in sv:
                    data_shifted[data_shifted == 0.0] = np.nan
                    data[data == 0.0] = np.nan

                pressure = dss.pressure
                pressure_interp = cf.interpolate_pressure(pressure)

                if np.sum(~np.isnan(data_shifted)) > 0:
                    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(14, 10), sharex=True, sharey=True)

                    # iterate through each profile and plot the profile lines
                    for pt in ptimes:
                        pt_idx = np.where(dss.profile_time.values == pt)[0]
                        non_nans = np.where(np.invert(np.isnan(data[pt_idx])))[0]
                        ax1.plot(data[pt_idx][non_nans], pressure_interp[pt_idx][non_nans], color='gray')  # plot lines

                        non_nans_shift = np.where(np.invert(np.isnan(data_shifted[pt_idx])))[0]
                        ax2.plot(data_shifted[pt_idx][non_nans_shift], pressure_interp[pt_idx][non_nans_shift], color='gray')  # plot lines

                    # add points
                    ax1.scatter(data, pressure_interp, color='gray', s=20, zorder=5)
                    ax2.scatter(data_shifted, pressure_interp, color='gray', s=20, zorder=5)

                    ax1.invert_yaxis()
                    ax1.set_ylabel('Pressure (dbar)')

                    ax1.ticklabel_format(useOffset=False)  # don't use scientific notation for ticks

                    try:
                        units = data.units
                    except AttributeError:
                        units = 'no_attributes'

                    ax1.set_xlabel(f'{data.name} ({units})')
                    ax2.set_xlabel(f'{data_shifted.name} ({units})')

                    opt_shift = pd.unique(dss[f'{sv.split("_shifted")[0]}_optimal_shift'].values).tolist()
                    for i, x in enumerate(opt_shift):
                        if ~np.isnan(x):
                            opt_shift[i] = int(x)

                    ax1.set_title('No shift')
                    ax2.set_title(f'Shifted: {opt_shift} seconds')
                    ttl = f'{deploy} {t0str} to {t1str}'
                    fig.suptitle(ttl)

                    sfile = os.path.join(savedir, save_filename)
                    plt.savefig(sfile, dpi=300)
                    plt.close()


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/gliders/qartod_qc/from_erddap/ru30-20210716T1804-profile-sci-rt.nc'
    save_directory = '/Users/garzio/Documents/rucool/gliders/qartod_qc/from_erddap/plotting'
    profile_group_n = 10
    main(ncfile, save_directory, profile_group_n)
