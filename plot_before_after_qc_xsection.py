#!/usr/bin/env python

"""
Author: lgarzio on 2/8/2022
Last modified: lgarzio on 2/16/2022
Plot cross-sections for glider science data defined by the user, before and after QC flags are applied.
Data are plotted from sci-profile datasets downloaded from RUCOOL's glider ERDDAP server using download_dataset.py
"""

import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable
import functions.common as cf
plt.rcParams.update({'font.size': 14})


def plot_xsection(figure, axis, x, y, z, cmap, title):
    xc = axis.scatter(x, y, c=z, cmap=cmap, s=12)

    # format colorbar
    divider = make_axes_locatable(axis)
    cax = divider.new_horizontal(size='5%', pad=0.1, axes_class=plt.Axes)
    figure.add_axes(cax)
    try:
        units = z.units
    except AttributeError:
        units = 'no_attributes'
    cb = plt.colorbar(xc, cax=cax, label=f'{z.name} ({units})')

    # add y-axis label
    axis.set_ylabel('Depth (m)')

    # format x-axis
    xfmt = mdates.DateFormatter('%Y\n%m-%d')
    axis.xaxis.set_major_formatter(xfmt)
    axis.xaxis.set_tick_params()

    axis.set_title(title, fontsize=16)


def flatten(lst):
    return [item for sublist in lst for item in sublist]


def main(ncf, sdir, inst_list, dr):
    fname = ncf.split('/')[-1].split('.nc')[0]
    deploy = f'{fname.split("-")[0]}-{fname.split("-")[1]}'

    ds = xr.open_dataset(ncf)
    ds = ds.swap_dims({'obs': 'time'})
    ds = ds.sortby(ds.time)

    if dr:
        ds = ds.sel(time=slice(dr[0], dr[1]))

    t0save = pd.to_datetime(np.nanmin(ds.time.values)).strftime('%Y%m%dT%H%M')
    t1save = pd.to_datetime(np.nanmax(ds.time.values)).strftime('%Y%m%dT%H%M')
    depth = ds.depth

    savedir = os.path.join(sdir, deploy, fname, 'xsection_before_after_qc')
    os.makedirs(savedir, exist_ok=True)

    varlist = []
    for il in inst_list:
        varlist.append(cf.define_instrument_variables(il))

    varlist = flatten(varlist)

    cmaps = cf.colormaps()

    for cv in varlist:
        save_filename = f'{cv}_xsection_before_after_qc_{t0save}-{t1save}.png'
        fig, (ax1, ax2) = plt.subplots(2, figsize=(14, 12), sharex=True, sharey=True)

        data = ds[cv]

        # in some cases, ERDDAP doesn't set the metadata/fill values, so get rid of any possible fill values
        data[data > 10000] = np.nan

        data_count = int(np.sum(~np.isnan(data)))

        # plot data without QC applied
        plot_xsection(fig, ax1, ds.time, depth, data, cmaps[cv], 'without QC')

        # find the qc summary variables and hysteresis test
        qc_vars = [x for x in ds.data_vars if f'{cv}_qartod_summary_flag' in x]
        if cv in ['conductivity', 'temperature']:
            qc_vars.append(f'{cv}_hysteresis_test')
        if cv in ['salinity', 'density']:
            qc_vars.append('conductivity_hysteresis_test')
            qc_vars.append('temperature_hysteresis_test')

        # apply QC variables to dataset and plot again
        qc_count = 0
        for qv in qc_vars:
            qc_idx = np.where(np.logical_or(ds[qv].values == 3, ds[qv].values == 4))[0]
            if len(qc_idx) > 0:
                data[qc_idx] = np.nan
                qc_count += len(qc_idx)

        percent = np.round(qc_count / data_count * 100)
        plot_xsection(fig, ax2, ds.time, depth, data, cmaps[cv], f'with QC (removed {qc_count} ({int(percent)}%) data points)')

        ax1.invert_yaxis()
        ttl = f'Before and After QC: {deploy}'
        fig.suptitle(ttl, y=0.95)

        sfile = os.path.join(savedir, save_filename)
        plt.savefig(sfile, dpi=300)
        plt.close()


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/gliders/qartod_qc/from_erddap/maracoos_02-20211020T1322-profile-sci-rt.nc'
    save_directory = '/Users/garzio/Documents/rucool/gliders/qartod_qc/from_erddap/plotting'
    instrument_list = ['ctd', 'do']  # ['ctd', 'do']
    date_range = None  # [dt.datetime(2021, 2, 28, 0, 0), dt.datetime(2021, 3, 5, 0, 0)]  # None
    main(ncfile, save_directory, instrument_list, date_range)
