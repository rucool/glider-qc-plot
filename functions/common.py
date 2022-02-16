#! /usr/bin/env python3

"""
Author: Lori Garzio on 1/26/2022
Last modified: 1/26/2022
"""
from erddapy import ERDDAP
import cmocean as cmo


def colormaps():
    cmaps = dict(conductivity='jet',
                 temperature=cmo.cm.thermal,
                 salinity=cmo.cm.haline,
                 density=cmo.cm.dense,
                 oxygen_concentration=cmo.cm.oxy,
                 oxygen_saturation='viridis')
    return cmaps


def define_instrument_variables(instrument):
    variables = dict(ctd=['conductivity', 'temperature', 'salinity', 'density'],
                     do=['oxygen_concentration', 'oxygen_saturation']
                     )
    return variables[instrument]


def flag_defs():
    flag_definitions = dict(unknown=dict(value=2, color='cyan'),
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
