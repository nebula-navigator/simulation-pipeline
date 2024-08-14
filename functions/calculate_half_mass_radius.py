#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:55:35 2024

@author: sohaib
"""

def calculate_half_mass_radius(data, mass_col):
    if mass_col not in data.columns or 'r' not in data.columns:
        raise ValueError(f"Data must contain '{mass_col}' and 'r' columns.")
    
    # Sort data by radial distance
    sorted_data = data.sort_values(by='r')
    print("sorted data half mass radius:", sorted_data)
    # Calculate the total mass
    total_mass = sorted_data[mass_col].sum()
    
    # Calculate cumulative mass
    sorted_data['cumulative_mass'] = sorted_data[mass_col].cumsum()
    print("total mass:",total_mass)
    # Determine radius where cumulative mass reaches half of the total mass
    half_mass = total_mass / 2.
    half_mass_radius = sorted_data.loc[sorted_data['cumulative_mass'] >= half_mass, 'r'].iloc[0]
    
    return half_mass_radius