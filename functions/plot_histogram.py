#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:59:21 2024

@author: sohaib
"""
import matplotlib.pyplot as plt
import seaborn as sns
import config




def plot_histogram(column_name):
    data=config.data
    if column_name not in data.columns:
        print(f"Error: Column '{column_name}' not found in data.")
        return
    plt.figure(figsize=(10, 6))
    sns.histplot(data[column_name].dropna(), kde=True)
    plt.title(f'Histogram of {column_name}')
    plt.xlabel(column_name)
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()