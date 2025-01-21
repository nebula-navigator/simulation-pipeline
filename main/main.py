#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:03:46 2024

@author: sohaib
"""
import sys
import os
import pandas as pd
import numpy as np
import config


from functions import (
    plot_histogram,
    plot_scatter,
    plot_scatter_3d,
    plot_boxplot,
    plot_pairplot,
    plot_violin,
    plot_heatmap,
    plot_3d_positions,
    plot_event_frequencies,
    plot_cumulative_number_vs_stellar_types,
    plot_position_vs_stellar_type,
    get_luminosity,
    plot_cumulative_counts_vs_r,
    plot_distribution_of_masses,
    calculate_core_radius,
    calculate_half_mass_radius,
    decode_history,
    get_cumulative_counts, history
)

from functions.suggest_visualizations import suggest_visualizations
from functions.history import history
from functions.historyplots import historyplots

# Global dictionary to store dataframes
data_dict = {}
current_data = None
column_names = [
    "im", "r", "vr", "vt", "u", "ikb", "a", "e", "ik1", "ik2", "sm1", "sm2", 
    "popId1", "popId2", "timenr", "idd1", "idd2", "lum1", "lum2", "rad1", 
    "rad2", "hist1", "hist2", "spin1", "spin2", "mu1", "mu2", "mv1", "mv2", 
    "mb1", "mb2", "mi1", "mi2", "mtr1", "mtr2", "star_type", "binary_star"
]

# k value dictionary
config.k_values = {
    0: 'Low-mass main-sequence star',
    1: 'Main-sequence star',
    2: 'Hertzsprung-gap star',
    3: 'First giant branch star',
    4: 'Core helium burning star',
    5: 'Early asymptotic giant branch star',
    6: 'Thermally pulsing asymptotic giant branch star',
    7: 'Naked helium star MS',
    8: 'Naked helium star Hertzsprung gap',
    9: 'Naked helium star giant branch',
    10: 'Helium white dwarf',
    11: 'Carbon-oxygen white dwarf',
    12: 'Oxygen-neon white dwarf',
    13: 'Neutron star',
    14: 'Black hole',
    15: 'Massless supernova'
}

k_values = config.k_values

def classify_stars(row):
    """Classify stars based on ik1 and ik2 values."""
    star_types = []
    if row['ikb'] == 0:  # Single star
        if pd.notna(row['ik1']):
            star_types.append(k_values.get(row['ik1'], 'Unknown'))
    else:  # Binary star
        if pd.notna(row['ik1']):
            star_types.append(k_values.get(row['ik1'], 'Unknown'))
        if pd.notna(row['ik2']):
            star_types.append(k_values.get(row['ik2'], 'Unknown'))
    return ', '.join(star_types) if star_types else 'Unknown'

def load_data(file_path):
    """Load data from a file and store it in the global dictionary."""
    file_key = os.path.basename(file_path) 
    data = pd.read_csv(file_path, names=column_names, delim_whitespace=True)
    config.data = data

    
    # Calculating projected positions if 'r' column is available
    if 'r' in data.columns:
        r = data['r'].values
        rpx, rpy, rpz = [], [], []
        
        for radius in r:
            cos_alpha_projection = np.random.uniform(-1., 1.)
            beta_projection = np.random.uniform(0., 1.) * 2 * np.pi
            gamma_projection = np.random.uniform(0., 1.) * 2 * np.pi
            
            alpha_projection = np.arccos(cos_alpha_projection)
            cos_beta_projection = np.cos(beta_projection)
            sin_beta_projection = np.sin(beta_projection)
            sin_alpha_projection = np.sin(alpha_projection)
            
            rpx_val = radius * cos_beta_projection * sin_alpha_projection
            rpy_val = radius * sin_beta_projection * sin_alpha_projection
            rpz_val = radius * cos_alpha_projection
            
            if abs(rpz_val) < 10**-15:
                rpz_val = 0.0
            if abs(rpx_val) < 10**-15:
                rpx_val = 0.0
            if abs(rpy_val) < 10**-15:
                rpy_val = 0.0
                
            rpx.append(rpx_val)
            rpy.append(rpy_val)
            rpz.append(rpz_val)
        
   
        data[f'posx'] = rpx
        data[f'posy'] = rpy
        data[f'posz'] = rpz
        
    
    
    # star_type column based on ik1 and ik2
    data['star_type'] = data.apply(classify_stars, axis=1)
    
    include = 'both'  
    cumulative_df = get_cumulative_counts(include)
    
   
    print("\nCumulative Counts DataFrame:")
    print(cumulative_df)
    
    list_columns(data)
    
    # Storing the dataframe in the global dictionary
    data_dict[file_key] = data
    
    return file_key

def list_data_files():
    """List all data files currently loaded."""
    if not data_dict:
        print("\nNo data files loaded.")
    else:
        print("\nLoaded data files:")
        for index, file_name in enumerate(data_dict.keys(), start=1):
            print(f"{index}. {file_name}")

def select_data_file():
    """Prompt the user to select a data file from the loaded ones."""
    list_data_files()
    if data_dict:
        try:
            choice = int(input("\nSelect a data file by number: "))
            file_name = list(data_dict.keys())[choice - 1]
            print(f"Selected file: {file_name}")
            global current_data
            current_data = data_dict[file_name]  # Update current_data based on user selection
            return file_name
        except (IndexError, ValueError):
            print("Invalid choice. Please enter a valid number.")
            return None
    return None


def list_columns(data):
    """List columns for the currently selected data."""
    print("\nAvailable columns for plotting:")
    for index, name in enumerate(data.columns, start=1):
        print(f"{index}. {name}")


def main_menu():
    global current_data
    exit_program = False  

    while not exit_program:
        print("\nChoose Analysis Type:")
        print("1. History File Analysis")
        print("2. Snap File Analysis")
        print("3. History Plot Generator")
        print("4. Exit")
        
        analysis_choice = input("Enter your choice (1, 2, or 3): ")
        
        if analysis_choice == '1':
            print("Opening history file analysis...")
            history()  
            print("History file analysis complete. Redirecting to analysis menu...")
            continue  

        elif analysis_choice == '2':
            print("Proceeding to snap file analysis menu.")
            while True:
                print("\nMain Menu:")
                print("1. Load Data Files")
                print("2. List Loaded Data Files")
                print("3. Select Data File for Operations")
                print("4. Suggest Visualization Types")
                print("5. Return to Analysis Menu")
                
                choice = input("Enter your choice (1-5): ")
                
                if choice == '1':
                    file_input = input("Enter the paths to data files, separated by commas: ")
                    file_paths = [path.strip() for path in file_input.split(',')]
                    for file_path in file_paths:
                        if os.path.exists(file_path):
                            file_key = load_data(file_path)
                            print(f"Data loaded successfully from: {file_path}")
                        else:
                            print(f"File not found: {file_path}")
                    print("\nData loading complete.")
                
                    if data_dict:
                        print("\nPlease select a data file to work with:")
                        select_data_file()
                        
                elif choice == '2':
                    list_data_files()
                elif choice == '3':
                    file_name = select_data_file()
                    if file_name:
                        data = data_dict[file_name]  
                        list_columns(current_data)
                elif choice == '4':
                    if current_data is not None:
                        suggest_visualizations(current_data)
                    else:
                        print("No data file selected. Please load and select a data file first.")
                elif choice == '5':
                    print("Returning to Analysis Menu.")
                    break 
                else:
                    print("Invalid choice. Please enter a number between 1 and 5.")
        elif analysis_choice == '3':
            historyplots()
            print("Plots Generated...Redirecting")
            continue
        elif analysis_choice == '4':
            print("Exiting the program.")
            exit_program = True 
        else:
            print("Invalid choice.")

    print("Come back soon. Let me know when your IMBH is ripe enough to be analyzed!")




if __name__ == "__main__":
    main_menu()
