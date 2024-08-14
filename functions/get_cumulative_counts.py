#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:54:37 2024

@author: sohaib
"""
import pandas as pd 
import config



def get_cumulative_counts(include='both'):
    
    data=config.data
    k_values=config.k_values
    cumulative_counts = {stype: 0 for stype in k_values.values()}

    # Initializing counts for single and binary
    single_counts = {stype: 0 for stype in k_values.values()}
    binary_counts = {stype: 0 for stype in k_values.values()}

    # Single stars
    if include in ['single', 'both']:
        single_data = data[data['ikb'] == 0]
        single_stars = single_data['ik1'].map(k_values).value_counts()
        
        for stype, count in single_stars.items():
            single_counts[stype] += count

    # Binary stars
    if include in ['binary', 'both']:
        binary_data = data[data['ikb'] != 0]
        binary_stars_1 = binary_data['ik1'].map(k_values).value_counts()
        binary_stars_2 = binary_data['ik2'].map(k_values).value_counts()
        
        for stype, count in binary_stars_1.items():
            binary_counts[stype] += count
        for stype, count in binary_stars_2.items():
            binary_counts[stype] += count

    # Combined counts for 'both'
    if include == 'both':
        for stype in cumulative_counts.keys():
            cumulative_counts[stype] = single_counts[stype] + binary_counts[stype]
    elif include == 'single':
        cumulative_counts.update(single_counts)
    elif include == 'binary':
        cumulative_counts.update(binary_counts)

    # Converting cumulative counts to a DataFrame
    cumulative_df = pd.DataFrame(list(cumulative_counts.items()), columns=['Stellar Type', 'Cumulative Count'])
    return cumulative_df