#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:02:15 2024

@author: sohaib
"""
import matplotlib.pyplot as plt
import seaborn as sns

def plot_heatmap(data):
   
    numeric_data = data.select_dtypes(include=[float, int])

    # Computing the correlation matrix
    corr_matrix = numeric_data.corr()

  
    plt.figure(figsize=(14, 10))

 
    sns.heatmap(
        corr_matrix, 
        annot=True, 
        fmt=".2f",  
        cmap='coolwarm', 
        cbar_kws={'shrink': .8},  
        linewidths=0.5,  
        square=True,  
        annot_kws={"size": 10}  
    )

    
    plt.title('Heatmap of Correlation Matrix', fontsize=16)

  
    plt.show()
