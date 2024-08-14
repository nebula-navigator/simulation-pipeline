#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:56:55 2024

@author: sohaib
"""

L_sun = 3.828e26  # Solar luminosity in Watts
k = 0.1
n = 0.5

def calculate_core_radius(luminosity, k=k, n=n):
    return k * (luminosity / L_sun) ** (1 / n)