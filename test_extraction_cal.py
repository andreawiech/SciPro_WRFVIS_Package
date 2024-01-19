#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 15:39:12 2024

@author: matildaachaab
"""

import numpy as np
import extration_and_calculation as ec



time_index = 10
lon =11
lat = 45

#calling functions to be tested
df = ec.extration_skewT_variables(time_index, lon, lat)
temp, pressure, mixing_ratio, geo_hght= ec.convert_var_to_actual_values(df)
dewpoint = ec.calculation_dewpoint(temp,pressure,mixing_ratio)


def test_extraction():
    df = ec.extration_skewT_variables(time_index, lon, lat)
    # Assuming variable_name and variable_units are correct attributes
    var_name = ['P', 'T', 'QVAPOR', 'U', 'V', 'PB', 'PHB', 'PH']
    var_unit = ['Pa', 'Unknown', 'kg kg-1', 'm s-1', 'm s-1', 'Pa', 'm2 s-2', 'm2 s-2']

    for i, variable in enumerate(var_name):
        assert 'variable_name' in df[variable].attrs, f"Attribute 'variable_name' not found for {variable}"
        assert 'units' in df[variable].attrs, f"Attribute 'variable_units' not found for {variable}"

        assert df[variable].attrs['variable_name'] == variable, f"Incorrect variable name for {variable}"
        assert df[variable].attrs['units'] == var_unit[i], f"Incorrect variable unit for {variable}"
        
 
        

def test_dewpoint(df):
    
    assert dewpoint.units == 'degree_Celsius'
    assert np.all(dewpoint.magnitude < 20.0), "dewpoint is not strictly decreasing"
    


def test_parcel_prof_temp(temp,pressure,dewpoint):
    
     prof = ec.parcel_profie(temp,pressure,dewpoint) 
     
     # Check if each element is strictly smaller than its predecessor
     for i in range(1, len(prof)):
        assert prof[i] < prof[i-1], "prof temp is not strictly decreasing"
     


def test_levels_skewt(pressure,temp,dewpoint):
    lcl_pressure, lcl_temperature,lfc_pressure, lfc_temperature = \
        ec.attris_of_skewT(pressure,temp,dewpoint)
    
    # Assert that LCL and LFC have the same pressure units
    assert lcl_pressure.units == lfc_pressure.units

    # Assert that LCL and LFC have the same temperature units
    assert lcl_temperature.units == lfc_temperature.units
    
    

    
    

