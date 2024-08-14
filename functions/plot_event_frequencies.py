#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:02:51 2024

@author: sohaib
"""

def plot_event_frequencies(column_name):
    if column_name not in data.columns:
        print(f"Error: Column '{column_name}' not found in data.")
        return
    plt.figure(figsize=(10, 6))
    sns.countplot(x=data[column_name].dropna())
    plt.title(f'Event Frequencies in {column_name}')
    plt.xlabel(column_name)
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()