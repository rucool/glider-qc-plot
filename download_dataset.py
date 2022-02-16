#!/usr/bin/env python

"""
Author: Lori Garzio on 1/26/2022
Last modified: 2/16/2022
Download a user-specified netCDF dataset from RUCOOL's glider ERDDAP server and save to a local directory
"""

import os
import functions.common as cf


def main(dsid, sdir, test, varlist, tr):
    os.makedirs(sdir, exist_ok=True)
    if test:
        erddap_server = 'http://slocum-test.marine.rutgers.edu//erddap'
    else:
        erddap_server = 'http://slocum-data.marine.rutgers.edu//erddap'

    kwargs = dict()
    if varlist:
        kwargs['variables'] = varlist
    if tr:
        kwargs['constraints'] = tr
    ds = cf.get_erddap_dataset(erddap_server, dsid, **kwargs)
    fname = f'{dsid}.nc'
    ds.to_netcdf(os.path.join(sdir, fname))


if __name__ == '__main__':
    dataset_id = 'maracoos_02-20211020T1322-profile-sci-rt'
    savedir = '/Users/garzio/Documents/rucool/gliders/qartod_qc/from_erddap'
    test = True  # True False
    variable_list = None
    time_range = None
    main(dataset_id, savedir, test, variable_list, time_range)
