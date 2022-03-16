#!/usr/bin/env python

"""
Author: lgarzio on 3/16/2022
Last modified: lgarzio on 3/16/2022
Plot time series of the optimal time shifts
Data are plotted from sci-profile datasets downloaded from RUCOOL's glider ERDDAP server using download_dataset.py
"""

import numpy as np
import pandas as pd
import xarray as xr
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
plt.rcParams.update({'font.size': 12})


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

    savedir = os.path.join(sdir, deploy, fname, f'timeseries_shift')
    os.makedirs(savedir, exist_ok=True)

    shift_vars = [v for v in list(ds.data_vars) if '_optimal_shift' in v]

    for sv in shift_vars:
        save_filename = f'{sv}_{t0save}-{t1save}.png'
        fig, ax = plt.subplots(figsize=(12, 8))

        data = ds[sv]

        # in some cases, ERDDAP doesn't set the metadata/fill values, so get rid of any possible fill values
        data[data > 10000] = np.nan

        # add points
        ax.scatter(data.time.values, data.values, color='k', s=5)

        try:
            units = data.units
        except AttributeError:
            units = 'no_attributes'

        ax.set_ylabel(f'{data.name} ({units})')
        plt.ylim([0, 60])
        ttl = f'{deploy}: optimal time shift calculated by glider segment'
        ax.set_title(ttl)

        # format x-axis
        xfmt = mdates.DateFormatter('%Y\n%m-%d')
        ax.xaxis.set_major_formatter(xfmt)
        ax.xaxis.set_tick_params()

        sfile = os.path.join(savedir, save_filename)
        plt.savefig(sfile, dpi=300)
        plt.close()


if __name__ == '__main__':
    ncfile = '/Users/garzio/Documents/rucool/gliders/qartod_qc/from_erddap/ru30-20210716T1804-profile-sci-rt.nc'
    save_directory = '/Users/garzio/Documents/rucool/gliders/qartod_qc/from_erddap/plotting'
    date_range = None  # [dt.datetime(2021, 2, 28, 0, 0), dt.datetime(2021, 3, 5, 0, 0)]  # None
    main(ncfile, save_directory, date_range)
