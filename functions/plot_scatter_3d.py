#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:00:55 2024

@author: sohaib
"""

import matplotlib.pyplot as plt
import config



def plot_scatter_3d(x_column, y_column, z_column,current_data, color_column=None):
    data=current_data
    if x_column not in data.columns or y_column not in data.columns or z_column not in data.columns:
        print(f"Error: Columns '{x_column}', '{y_column}', or '{z_column}' not found in data.")
        return
    
    x = data[x_column].dropna().values
    y = data[y_column].dropna().values
    z = data[z_column].dropna().values
    
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    if color_column and color_column in data.columns:
        c = data[color_column].dropna().values
        scatter = ax.scatter(x, y, z, c=c, cmap='viridis')
        fig.colorbar(scatter, ax=ax, label=color_column)
    else:
        ax.scatter(x, y, z, color='blue')

    ax.set_title(f'3D Scatter Plot of {x_column}, {y_column}, and {z_column}')
    ax.set_xlabel(x_column)
    ax.set_ylabel(y_column)
    ax.set_zlabel(z_column)
    plt.grid(True)
    plt.show()