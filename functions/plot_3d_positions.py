#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:02:34 2024

@author: sohaib
"""
import numpy as np
import matplotlib.pyplot as plt
import config



def plot_3d_positions(current_data,color_column=None, color_min=None, color_max=None):
    data=current_data
    k_values=config.k_values
    if 'r' not in data.columns:
        print("Error: Column 'r' not found in data.")
        return
    
    filtered_data = data
    if color_column and color_min is not None and color_max is not None:
        filtered_data = data[(data[color_column] >= float(color_min)) & (data[color_column] <= float(color_max))]
    
    r = filtered_data['r'].values
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

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    if color_column and color_column in data.columns:
        c = filtered_data[color_column].values
        scatter = ax.scatter(rpx, rpy, rpz, c=c, cmap='viridis')
        fig.colorbar(scatter, ax=ax, label=color_column)
    else:
        ax.scatter(rpx, rpy, rpz, color='blue')

    ax.set_title('3D Projected Positions of Stars')
    ax.set_xlabel('X (parsecs)')
    ax.set_ylabel('Y (parsecs)')
    ax.set_zlabel('Z (parsecs)')
    plt.grid(True)
    plt.show()