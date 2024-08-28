import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import config
from PyAstronomy import pyasl
import math


def calculate_temperature(bv):
    
    if math.isnan(bv):
        return None  
    
    b = pyasl.BallesterosBV_T()
    teff = b.bv2T(bv)
    return teff


def plothr():
    
    data = config.data
    k_values = config.k_values

    
    available_types = sorted(k_values.keys())
    print("Available stellar types (k-values):")
    for i, stype in enumerate(available_types):
        print(f"{i}: {k_values[stype]}")  

    
    user_input = input("Enter the numbers of the stellar types you want to include, separated by commas (e.g., 0,1,3), or press Enter to include all: ")
    if user_input.strip() == '':
        selected_types = available_types
    else:
        selected_types = [int(x.strip()) for x in user_input.split(',') if x.strip().isdigit()]

    print(f"Selected stellar types: {[k_values.get(stype, 'Unknown') for stype in selected_types]}")

    # Handle k=0 cases
    if 0 in selected_types:
        data.loc[data['sm1'] == 0, 'ik1'] = np.nan
        data.loc[data['sm2'] == 0, 'ik2'] = np.nan

    # masks
    mask1 = data['sm1'] != 0
    mask2 = data['sm2'] != 0

    # B - V color index
    data['bv1'] = data['mb1'] - data['mv1']
    data['bv2'] = data['mb2'] - data['mv2']
    

    data.loc[data['bv1'] == 0.0, 'bv1'] = np.nan
    data.loc[data['bv2'] == 0.0, 'bv2'] = np.nan
   
    data['temp1'] = data['bv1'].apply(lambda bv: calculate_temperature(bv) if not pd.isna(bv) else np.nan)
    data['temp2'] = data['bv2'].apply(lambda bv: calculate_temperature(bv) if not pd.isna(bv) else np.nan)

    
    min_temp = min(data[['temp1', 'temp2']].min())
    max_temp = max(data[['temp1', 'temp2']].max())

    
    colors = plt.cm.tab20(np.linspace(0, 1, len(selected_types)))
    color_map = dict(zip(selected_types, colors))

   
    combined_bv = []
    combined_magnitude = []
    combined_temp = []
    labels = []
    marker_colors = []
    combined_index = []

    for stype in selected_types:
        # Apply masks for each type
        mask_stype1 = (data['ik1'] == stype)
        mask_stype2 = (data['ik2'] == stype)

        bv_stype1 = data.loc[mask1 & mask_stype1, 'bv1']
        mv_stype1 = data.loc[mask1 & mask_stype1, 'mv1']
        temp_stype1 = data.loc[mask1 & mask_stype1, 'temp1']

        bv_stype2 = data.loc[mask2 & mask_stype2, 'bv2']
        mv_stype2 = data.loc[mask2 & mask_stype2, 'mv2']
        temp_stype2 = data.loc[mask2 & mask_stype2, 'temp2']

        # Append data for plotting
        combined_bv.extend(bv_stype1.tolist())
        combined_bv.extend(bv_stype2.tolist())
        combined_magnitude.extend(mv_stype1.tolist())
        combined_magnitude.extend(mv_stype2.tolist())
        combined_temp.extend(temp_stype1.tolist())
        combined_temp.extend(temp_stype2.tolist())
        labels.extend([k_values.get(stype, 'Unknown')] * (len(bv_stype1) + len(bv_stype2)))
        marker_colors.extend([color_map[stype]] * (len(bv_stype1) + len(bv_stype2)))
        combined_index.extend(data.index[mask1 & (data['ik1'] == stype)].tolist() + data.index[mask2 & (data['ik2'] == stype)].tolist())

   
    combined_bv = pd.Series(combined_bv)
    combined_magnitude = pd.Series(combined_magnitude)
    combined_temp = pd.Series(combined_temp)
    combined_index = pd.Series(combined_index)

    
    marker_size = 25  

    fig, ax1 = plt.subplots(figsize=(12, 8))

    
    scatter_plot = ax1.scatter(combined_bv, combined_magnitude, 
                               c=marker_colors, alpha=0.7, 
                               s=marker_size, marker='.', edgecolors='none', picker=True)

    
    ax2 = ax1.twiny()
    ax2.set_xlim(ax1.get_xlim())
    
    
    num_ticks = len(ax1.get_xticks())  
    temperature_ticks = np.linspace(min_temp, max_temp, num_ticks)  

    
    ax2.set_xticks(np.linspace(min(combined_bv), max(combined_bv), num_ticks)) 
    ax2.set_xticklabels([f"{int(temp)} K" for temp in reversed(temperature_ticks)])  
    ax2.set_xlabel('Temperature (K)')

    # Add labels, title, and grid
    ax1.set_xlabel('B - V Color Index')
    ax1.set_ylabel('Absolute Visual Magnitude (Mv)')
    ax1.set_title('Hertzsprung-Russell Diagram')
    ax1.invert_yaxis() 
    ax1.grid(True)
    
    
    handles = [plt.Line2D([0], [0], marker='.', color='w', markerfacecolor=color_map[stype], markersize=10, label=k_values.get(stype, 'Unknown')) for stype in selected_types]
    ax1.legend(title='Stellar Type', handles=handles, loc='best')

    plt.tight_layout()
    plt.show()

    print("Updated data with calculated temperatures for components 'ik1' and 'ik2':")
    print(data[['bv1', 'bv2', 'temp1', 'temp2']].head())
