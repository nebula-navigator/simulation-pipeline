#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:58:45 2024

@author: sohaib
"""

def get_luminosity(mass1, mass2, lum1, lum2):
    if mass1 >= mass2:
        return lum1
    else:
        return lum2