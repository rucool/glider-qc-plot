#! /usr/bin/env python3

"""
Author: Lori Garzio on 1/26/2022
Last modified: 3/16/2022
"""
import xarray as xr
import numpy as np
from erddapy import ERDDAP
import cmocean as cmo


def colormaps():
    cmaps = dict(conductivity='jet',
                 temperature=cmo.cm.thermal,
                 salinity=cmo.cm.haline,
                 density=cmo.cm.dense,
                 oxygen_concentration=cmo.cm.oxy,
                 oxygen_saturation='viridis',
                 oxygen_concentration_shifted=cmo.cm.oxy,
                 oxygen_saturation_shifted='viridis',
                 sbe41n_ph_ref_voltage_shifted=cmo.cm.matter
                 )
    return cmaps


def define_instrument_variables(instrument):
    variables = dict(ctd=['conductivity', 'temperature', 'salinity', 'density'],
                     do=['oxygen_concentration', 'oxygen_saturation']
                     )
    return variables[instrument]


def flag_defs():
    flag_definitions = dict(not_evaluated=dict(value=2, color='cyan'),
                            suspect=dict(value=3, color='orange'),
                            fail=dict(value=4, color='red'))
    return flag_definitions


def flatten(lst):
    return [item for sublist in lst for item in sublist]


def get_erddap_dataset(server, ds_id, variables=None, constraints=None):
    e = ERDDAP(server=server,
               protocol='tabledap',
               response='nc')
    e.dataset_id = ds_id
    if constraints:
        e.constraints = constraints
    if variables:
        e.variables = variables
    ds = e.to_xarray()
    return ds


def interpolate_pressure(da, varname='pressure'):
    """
    Linear interpolate pressure.
    :param da: data in the form of an xarray dataset
    :param varname: variable name to interpolate, default is pressure
    :returns: interpolated pressure in the form of an xarray dataset
    """
    df = da.to_dataframe()
    if varname == 'depth':
        df = df.iloc[:, :-1]
    df[varname] = df[varname].interpolate(method='linear', limit_direction='both')
    interp_da = xr.DataArray(np.array(df[varname]).astype(da.dtype), coords=da.coords, dims=da.dims,
                             name=f'{varname}_interpolated', attrs=da.attrs.copy())
    return interp_da
