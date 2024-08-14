#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:03:11 2024

@author: sohaib
"""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import config



def plot_cumulative_number_vs_stellar_types(include):
    data=config.data
    k_values=config.k_values
    include = include.lower()
    valid_includes = ['single', 'binary', 'both']
    if include not in valid_includes:
        print(f"Error: Include must be one of {valid_includes}.")
        return

    cumulative_counts = {stype: 0 for stype in k_values.values()}

    # Initialize counts for single and binary
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

    # Converting cumulative counts to a DataFrame for plotting
    cumulative_counts_df = pd.DataFrame(list(cumulative_counts.items()), columns=['Stellar Type', 'Count'])
    
    # Plot the cumulative counts
    plt.figure(figsize=(12, 6))
    sns.barplot(x='Stellar Type', y='Count', data=cumulative_counts_df)
    plt.xticks(rotation=90)
    plt.title('Cumulative Number vs Stellar Types')
    plt.xlabel('Stellar Type')
    plt.ylabel('Cumulative Number')
    plt.grid(True)
    plt.show()