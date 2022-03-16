#!/usr/bin/env python

"""
Author: lgarzio on 3/16/2022
Last modified: lgarzio on 3/16/2022
Plot cross-sections for glider data, before and after time shifting
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


def main(ncf, sdir, dr):
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
    depth_interp = cf.interpolate_pressure(depth, varname='depth')

    savedir = os.path.join(sdir, deploy, fname, 'xsection_before_after_timeshift')
    os.makedirs(savedir, exist_ok=True)

    shifted_vars = [v for v in list(ds.data_vars) if '_shifted' in v]

    cmaps = cf.colormaps()

    for sv in shifted_vars:
        save_filename = f'{sv}_xsection_before_after_shift_{t0save}-{t1save}.png'
        fig, (ax1, ax2) = plt.subplots(2, figsize=(14, 12), sharex=True, sharey=True)

        data_shifted = ds[sv]
        data = ds[sv.split('_shifted')[0]]

        # in some cases, ERDDAP doesn't set the metadata/fill values, so get rid of any possible fill values
        data[data > 10000] = np.nan

        # convert 0.0 pH voltage values to nan
        if 'ph_ref_voltage' in sv:
            data_shifted[data_shifted == 0.0] = np.nan
            data[data == 0.0] = np.nan

        plot_xsection(fig, ax1, ds.time, depth_interp, data, cmaps[sv], 'No Shift')
        plot_xsection(fig, ax2, ds.time, depth_interp, data_shifted, cmaps[sv], 'Shifted')

        ax1.invert_yaxis()
        ttl = f'Before and After Time Shifting: {deploy}'
        fig.suptitle(ttl, y=0.95)

        sfile = os.path.join(savedir, save_filename)
        plt.savefig(sfile, dpi=300)
        plt.close()


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/gliders/qartod_qc/from_erddap/ru30-20210716T1804-profile-sci-rt.nc'
    save_directory = '/Users/garzio/Documents/rucool/gliders/qartod_qc/from_erddap/plotting'
    date_range = None  # [dt.datetime(2021, 2, 28, 0, 0), dt.datetime(2021, 3, 5, 0, 0)]  # None
    main(ncfile, save_directory, date_range)
