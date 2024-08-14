#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:00:00 2024

@author: sohaib
"""
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import config



def plot_scatter(x_column, y_column, continuous_column=None, color_min=None, color_max=None, categorical_column=None):
    data=config.data
    k_values=config.k_values
    if x_column not in data.columns or (y_column != 'cumulative_count' and y_column not in data.columns):
        print(f"Error: Columns '{x_column}' or '{y_column}' not found in data.")
        return
    
    # Always prompt for star type selection
    print("\nAvailable star types (based on k values):")
    for k, v in k_values.items():
        print(f"{k}: {v}")
    
    selected_ks = input("Enter the k values to include (comma-separated) or press Enter to include all: ")
    
    if selected_ks.strip():
        selected_ks = [int(k.strip()) for k in selected_ks.split(',')]
        # Filter data to include only selected star types
        filtered_data = data[data['ik1'].isin(selected_ks) | data['ik2'].isin(selected_ks)]
    else:
        filtered_data = data

    # Apply filter using the continuous variable column, if provided
    if continuous_column and continuous_column in filtered_data.columns:
        c = filtered_data[continuous_column]
        
        if color_min is not None and color_max is not None:
            mask = (c >= float(color_min)) & (c <= float(color_max))
            range_text = f'{continuous_column} range: [{color_min}, {color_max}]'
        elif color_min is not None:
            mask = (c >= float(color_min))
            range_text = f'{continuous_column} range: [{color_min}, max]'
        elif color_max is not None:
            mask = (c <= float(color_max))
            range_text = f'{continuous_column} range: [min, {color_max}]'
        else:
            mask = np.ones_like(c, dtype=bool)  # No filtering
            range_text = f'{continuous_column} range: [min, max]'
        
        filtered_data = filtered_data[mask]
    else:
        range_text = None  # No range filtering applied

    # Ask for optional limits on the categorical column
    if categorical_column and categorical_column in filtered_data.columns:
        sec_color_min = input(f"Enter minimum value for {categorical_column} (or press Enter to use default): ")
        sec_color_max = input(f"Enter maximum value for {categorical_column} (or press Enter to use default): ")
        sec_color_min = float(sec_color_min) if sec_color_min else None
        sec_color_max = float(sec_color_max) if sec_color_max else None

        if sec_color_min is not None or sec_color_max is not None:
            c = filtered_data[categorical_column]
            if sec_color_min is not None and sec_color_max is not None:
                sec_mask = (c >= sec_color_min) & (c <= sec_color_max)
            elif sec_color_min is not None:
                sec_mask = (c >= sec_color_min)
            elif sec_color_max is not None:
                sec_mask = (c <= sec_color_max)
            else:
                sec_mask = np.ones_like(c, dtype=bool)
            
            filtered_data = filtered_data[sec_mask]
            y_values = filtered_data[y_column]
        else:
            y_values = filtered_data[y_column]
    else:
        y_values = filtered_data[y_column]
        sec_color_min = sec_color_max = None  # No secondary color limits

    # Convert categorical columns to categorical data type
    if categorical_column in filtered_data.columns:
        filtered_data[categorical_column] = filtered_data[categorical_column].astype('category')
    
    plt.figure(figsize=(10, 6))
    
   
    
    
        # Scatter plot with the categorical column as hue
    sns.scatterplot(
            x=filtered_data[x_column], 
            y=y_values,
            hue=filtered_data[categorical_column],  # Color based on the categorical column
            style=filtered_data['star_type'],  # Use different markers for star types
            palette='viridis',
            alpha=1, s=100,
            legend='full'
        )
        
    # Add text box showing the selected range of the continuous variable in the top-left corner
    if range_text:
        plt.gca().text(0.02, 0.98, range_text, transform=plt.gca().transAxes, 
                       fontsize=10, verticalalignment='top', horizontalalignment='left',
                       bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='lightyellow'))
    
    # Move the legend outside of the plot and adjust the layout
    plt.legend(title='Legend', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.title(f'Scatter Plot of {x_column} vs {y_column}' if y_column != 'cumulative_count' else 'Line Plot of Cumulative Counts')
    plt.xlabel(x_column)
    plt.ylabel(y_column)
    plt.grid(True)
    
    plt.show(block=False)
    plt.pause(0.1)  # Allow the plot to display
    
    # Prompt user if they want to save the filtered data
    save_data_prompt = input("Would you like to save the filtered data? (y/n): ").strip().lower()
    
    if save_data_prompt == 'y':
        file_name = input("Enter the file name to save the data (without extension): ").strip() + '.dat'
        metadata = f"# Filtered by {continuous_column} range: [{color_min}, {color_max}]\n"
        metadata += f"# Selected k-values: {selected_ks}\n"
        metadata += f"# x_column: {x_column}, y_column: {y_column}\n"
        metadata += f"# categorical_column: {categorical_column} with limits [{sec_color_min}, {sec_color_max}]\n"
        
        # Align columns
        column_widths = [max(len(str(value)) for value in [col] + filtered_data[col].tolist()) + 2 for col in filtered_data.columns]
        with open(file_name, 'w') as f:
            f.write("# Metadata\n")
            f.write(metadata)
            f.write("# Data\n")
            
            # Write column headers with proper alignment
            for col, width in zip(filtered_data.columns, column_widths):
                f.write(f"{col.ljust(width)}")
            f.write('\n')
            
            # Write data rows with proper alignment
            for _, row in filtered_data.iterrows():
                for col, width in zip(filtered_data.columns, column_widths):
                    f.write(f"{str(row[col]).ljust(width)}")
                f.write('\n')
        
        print(f"Data saved to {file_name}")
    else:
        print("Data not saved.")