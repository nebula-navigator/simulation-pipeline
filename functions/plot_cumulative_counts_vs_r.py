#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:57:55 2024

@author: sohaib
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import config
from .calculate_core_radius import calculate_core_radius
from .calculate_half_mass_radius import calculate_half_mass_radius
from .get_luminosity import get_luminosity



def plot_cumulative_counts_vs_r():
    data=config.data
    k_values=config.k_values
    include_options = ['single', 'binary', 'both']
    print("Select the type of stars to include:")
    print("1. Single Stars")
    print("2. Binary Stars")
    print("3. Both")

    choice = input("Enter your choice (1-3): ")
    
    if choice == '1':
        include = 'single'
    elif choice == '2':
        include = 'binary'
    elif choice == '3':
        include = 'both'
    else:
        print("Invalid choice. Defaulting to 'both'.")
        include = 'both'
    
    if include == 'single':
        filtered_data = data[data['ikb'] == 0]
    elif include == 'binary':
        filtered_data = data[data['ikb'] != 0]
    else:
        filtered_data = data
    
    print("Available stellar types:")
    for k, stype in k_values.items():
        print(f"{k}: {stype}")

    valid_types_input = input("Enter stellar type numbers to include, separated by commas (e.g., '1,2,3'): ")
    
    if valid_types_input:
        valid_types = [int(t.strip()) for t in valid_types_input.split(',')]
    else:
        valid_types = list(k_values.keys())

    filtered_data = filtered_data[
        (filtered_data['ik1'].isin(valid_types)) | 
        (filtered_data['ik2'].isin(valid_types))
    ]
    
    if include == 'single':
        half_mass_radius = calculate_half_mass_radius(filtered_data, 'sm1')
    else:
        filtered_data = filtered_data.copy()
        filtered_data['combined_mass'] = filtered_data['sm1'] + filtered_data['sm2']
        half_mass_radius = calculate_half_mass_radius(filtered_data, 'combined_mass')

    # Calculate core radius for binaries
    filtered_data['effective_lum'] = filtered_data.apply(
        lambda row: get_luminosity(row['sm1'], row['sm2'], row['lum1'], row['lum2']),
        axis=1
    )
    filtered_data['core_radius'] = calculate_core_radius(filtered_data['effective_lum'])
    
    core_radius_mean = filtered_data['core_radius'].mean()
    
    print(f"Half-Mass Radius: {half_mass_radius:.2f} pc")
    print(f"Mean Core Radius: {core_radius_mean:.2f} pc")
    
    min_radial = input("Enter minimum radial position (pc) or press Enter to skip: ")
    max_radial = input("Enter maximum radial position (pc) or press Enter to skip: ")
    
    if min_radial:
        min_radial = float(min_radial)
        filtered_data = filtered_data[filtered_data['r'] >= min_radial]
    
    if max_radial:
        max_radial = float(max_radial)
        filtered_data = filtered_data[filtered_data['r'] <= max_radial]
    
    if filtered_data.empty:
        print("No data available for the selected filters.")
        return
    
    cumulative_counts_df = pd.DataFrame()
    grouped = filtered_data.groupby(['r', 'ik1']).size().unstack(fill_value=0)
    
    print("Grouped DataFrame columns (stellar types):", grouped.columns)
    
    for stype, stype_name in k_values.items():
        if stype in valid_types:
            if stype in grouped.columns:
                cumulative_counts = grouped[stype].cumsum()

                # Ensure we do not drop to zero prematurely
                max_count = cumulative_counts.max()
                if max_count > 0:
                    constant_index = cumulative_counts[cumulative_counts == max_count].index[0]
                    cumulative_counts.loc[constant_index+1:] = 0  # Drop to zero after it becomes constant
                
                cumulative_counts_df[stype_name] = cumulative_counts
            else:
                print(f"Stellar type {stype_name} not present in data.")
                cumulative_counts_df[stype_name] = pd.Series(index=grouped.index, dtype=float).fillna(0).cumsum()
    
    r_values = cumulative_counts_df.index
    
    plt.figure(figsize=(10, 6))
    
    for stype_name in cumulative_counts_df.columns:
        plt.plot(r_values, cumulative_counts_df[stype_name], label=stype_name)
        
        # Add star type label at the end of each line
        end_x = r_values[-1]
        end_y = cumulative_counts_df[stype_name].iloc[-1]
        plt.text(end_x, end_y, stype_name, fontsize=8, verticalalignment='bottom', horizontalalignment='right')
    
    plt.axvline(x=half_mass_radius, color='r', linestyle='--', label=f'Half-Mass Radius: {half_mass_radius:.2f} pc')
    plt.axvline(x=core_radius_mean, color='b', linestyle='--', label=f'Mean Core Radius: {core_radius_mean:.2f} pc')
    
    plt.xlabel('Radial Position (pc)')
    plt.ylabel('Cumulative Count')
    plt.title('Cumulative Counts vs. Radial Position')
    plt.grid(True)
    plt.legend(title='Stellar Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.show()