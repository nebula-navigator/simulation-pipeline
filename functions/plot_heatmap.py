#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:02:15 2024

@author: sohaib
"""
import matplotlib.pyplot as plt
import seaborn as sns

def plot_heatmap(data):
    # Exclude non-numeric columns
    numeric_data = data.select_dtypes(include=[float, int])

    # Compute the correlation matrix
    corr_matrix = numeric_data.corr()

    # Set up the matplotlib figure
    plt.figure(figsize=(14, 10))

    # Create the heatmap with more readability features
    sns.heatmap(
        corr_matrix, 
        annot=True, 
        fmt=".2f",  # Limit the annotation to 2 decimal places
        cmap='coolwarm', 
        cbar_kws={'shrink': .8},  # Shrink the color bar to fit better
        linewidths=0.5,  # Add lines between cells for clarity
        square=True,  # Make the cells square-shaped
        annot_kws={"size": 10}  # Control the size of the annotations
    )

    # Add a title to the heatmap
    plt.title('Heatmap of Correlation Matrix', fontsize=16)

    # Show the plot
    plt.show()
