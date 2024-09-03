#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:03:30 2024

@author: sohaib
"""
import matplotlib.pyplot as plt
import seaborn as sns
import config



def plot_position_vs_stellar_type(current_data):
    data=current_data
    k_values=config.k_values
    include = input("Include 'single', 'binary', or 'both' stars: ").strip().lower()
    valid_includes = ['single', 'binary', 'both']
    
    if include not in valid_includes:
        print(f"Error: Include must be one of {valid_includes}.")
        return

    if include == 'single':
        valid_data = data[data['ikb'] == 0]
    elif include == 'binary':
        valid_data = data[data['ikb'] != 0]
    else:
        valid_data = data

    valid_data['Stellar Type'] = valid_data.apply(lambda row: k_values.get(row['ik1'], 'Unknown'), axis=1)
    
    plt.figure(figsize=(12, 6))
    sns.scatterplot(x='r', y='Stellar Type', data=valid_data, alpha=0.6)
    plt.title('Position (r) vs Stellar Type')
    plt.xlabel('Position (r)')
    plt.ylabel('Stellar Type')
    plt.grid(True)
    plt.show()