#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:51:32 2024

@author: sohaib
"""

def decode_history(hist_value):
    events = [
        "dynamical collision (inter-coll)",
        "binary merger (stellar evolution)",
        "mass transfer (increased the mass)",
        "interaction (bin-bin, or bin-sin)",
        "exchange in interaction",
        "went to binary (inter-binform)",
        "went to single (dissolution in stellar evolution)",
        "went to single (dissolution in interaction)",
        "separation changed by 1% (in interaction)",
        "separation changed by 10% (in interaction)",
        "eccentricity changed by 1% (in interaction)",
        "eccentricity changed by 10% (in interaction)",
        "dissolution in bin-sin interaction",
        "dissolution in bin-bin interaction",
        "dissolution in binary evolution",
        "exchange in bin-sin interaction",
        "exchange in bin-bin interaction",
        "dynamical collision in bin-sin",
        "dynamical collision in bin-bin (or hierarchical)"
    ]
    
    
    binary_string = format(hist_value, '020b')[::-1]  # 20-bit binary string, reversed
    
    # Decode events
    decoded_events = [
        events[i] for i, bit in enumerate(binary_string) if bit == '1' and i < len(events)
    ]
    return decoded_events