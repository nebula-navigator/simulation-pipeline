#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:01:25 2024

@author: sohaib
"""
import matplotlib.pyplot as plt
import seaborn as sns
import config




def plot_boxplot(column_name):
    data=config.data
    if column_name not in data.columns:
        print(f"Error: Column '{column_name}' not found in data.")
        return
    plt.figure(figsize=(10, 6))
    sns.boxplot(y=data[column_name].dropna())
    plt.title(f'Box Plot of {column_name}')
    plt.xlabel(column_name)
    plt.grid(True)
    plt.show()