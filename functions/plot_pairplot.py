#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:01:43 2024

@author: sohaib
"""
import matplotlib.pyplot as plt
import seaborn as sns
import config



def plot_pairplot(columns_list):
    data=config.data
    k_values=config.k_values
    if any(col not in data.columns for col in columns_list):
        print("Error: One or more columns not found in data.")
        return
    sns.pairplot(data[columns_list].dropna())
    plt.title('Pair Plot')
    plt.grid(True)
    plt.show()
