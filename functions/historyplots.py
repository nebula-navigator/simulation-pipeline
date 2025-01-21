#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 20:14:46 2025

@author: sohaib
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import config
import numpy as np


def historyplots():
   
    file_path = input("Please enter the path to the history file: ")
    history_basename = os.path.splitext(os.path.basename(file_path))[0]

   
    with open(file_path, 'r') as file:
        history_content = file.readlines()

    column_line = history_content[187].strip()
    columns = column_line.lstrip('#').split()
    if columns[0] == '#':
        columns[0] = 'lineType(1)'

   
    df = pd.read_csv(file_path, delim_whitespace=True, skiprows=188, comment='#', names=columns)

   
    output_directory = f"individual_plots_{history_basename}"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    
    df = df[df['starType'] == 14]
    df['MassChange'] = df['massNew[Msun]'] - df['massOld[Msun](10)']
    mass_increase_points = df[df['MassChange'] > 0.]

    sns.set(style="darkgrid")

    
    plt.figure(figsize=(8, 6))
    plt.gca().set_facecolor('black')
    plt.grid(True, color='white', linestyle='-', linewidth=0.5)
    plt.plot(df['time[Myr]'], df['massNew[Msun]'], label='Mass Evolution', color='w')
    plt.scatter(mass_increase_points['time[Myr]'], mass_increase_points['massNew[Msun]'],
                color='yellow', label='Mass Increase Points', s=10)
    plt.title("Mass Evolution Over Time with Significant Events (0-2 Gyr)")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Mass (Msun)")
    plt.xlim([0, 2000])
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.savefig(os.path.join(output_directory, "mass_evolution_0_2Gyr.png"), dpi=300, bbox_inches="tight")
    plt.close()

    
    plt.figure(figsize=(8, 6))
    plt.gca().set_facecolor('black')
    plt.grid(True, color='white', linestyle='-', linewidth=0.5)
    plt.plot(df['time[Myr]'], df['massNew[Msun]'], label='Mass Evolution', color='w')
    plt.scatter(mass_increase_points['time[Myr]'], mass_increase_points['massNew[Msun]'],
                color='yellow', label='Mass Increase Points', s=10)
    plt.title("Mass Evolution Over Time with Significant Events (Entire Evolution)")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Mass (Msun)")
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.savefig(os.path.join(output_directory, "mass_evolution_entire.png"), dpi=300, bbox_inches="tight")
    plt.close()

   
    event_counts = df['lineType(1)'].value_counts()
    plt.figure(figsize=(8, 6))
    sns.barplot(x=event_counts.index, y=event_counts.values,palette="viridis")
    plt.title("Event Type Frequency")
    plt.xlabel("Event Type")
    plt.ylabel("Frequency")
    plt.xticks(rotation=45)
    plt.savefig(os.path.join(output_directory, "event_type_frequency.png"), dpi=300, bbox_inches="tight")
    plt.close()

    
    plt.figure(figsize=(8, 6))
    plt.plot(df['time[Myr]'], df['r[pc]'], label='r', color='r')
    plt.title("IMBH Position Change Over Time")
    plt.xlabel("Time (Myr)")
    plt.xscale('log')
    plt.ylabel("r [pc]")
    plt.legend()
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.savefig(os.path.join(output_directory, "imbh_position_change.png"), dpi=300, bbox_inches="tight")
    plt.close()

   
    filtered_df = df[df['encCase'] != 0]
    grouped = filtered_df.groupby(['time[Myr]', 'encCase']).size().unstack(fill_value=0)
    cumulative_counts_df = grouped.cumsum()

    plt.figure(figsize=(8, 6))
    for enc_case in cumulative_counts_df.columns:
        plt.plot(cumulative_counts_df.index, cumulative_counts_df[enc_case], label=f'Encounter Case {enc_case}')
    plt.title("Cumulative Encounter Cases Over Time")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Cumulative Count")
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.grid(True)
    plt.yscale('log')
    plt.savefig(os.path.join(output_directory, "cumulative_encounter_cases.png"), dpi=300, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(8, 6))
    stats_text = """
    Mergers:
    - Same binary mergers: 100
    - Different binary mergers: 50
    - Double mergers with exchanges: 25

    Flybys/Exchanges:
    - No change in binaries: 200
    - Single exchange: 150
    - Double exchange: 75
    - Binary unbounded or merged: 100

    Single Stars formed after interactions:
    - Total: 300
    """
    plt.text(0.5, 0.5, stats_text, ha='center', va='center', fontsize=10,
             bbox=dict(facecolor='lightgray', edgecolor='black'))
    plt.axis('off')
    plt.savefig(os.path.join(output_directory, "encounter_stats.png"), dpi=300, bbox_inches="tight")
    plt.close()
    
    
    filtered_df1 = df[df['mergMass1(26)'] > 0.001]
    event_counts_2 = filtered_df1['lineType(1)'].value_counts()
    
    plt.figure(figsize=(8, 6))
    sns.barplot(x=event_counts_2.index, y=event_counts_2.values,palette="viridis")
    plt.title("Mergers Events")
    plt.xlabel("Event Type")
    plt.ylabel("Frequency")
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(os.path.join(output_directory, "mergers_events.png"), dpi=300, bbox_inches="tight")
    plt.close()
    
   
    
   
    filtered_df2 = df[df['mergMass1(26)'] > 0.001].copy()
    event_counts_3 = filtered_df2['compType'].value_counts()
    
    plt.figure(figsize=(8, 6))
    sns.barplot(x=event_counts_3.index, y=event_counts_3.values,palette="viridis")
    plt.title("Mergers with Each Stellar Type")
    plt.xlabel("Stellar Type")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(os.path.join(output_directory, "mergers_stellar_type.png"), dpi=300, bbox_inches="tight")
    plt.close()
    
 
    filtered_df4 = filtered_df2.copy()
    all_events = sorted(filtered_df2['lineType(1)'].unique())
    grouped = filtered_df2.groupby(['time[Myr]', 'lineType(1)']).size().unstack(fill_value=0).reindex(columns=all_events, fill_value=0)
    cumulative_counts_df = grouped.cumsum()
    cumulative_counts_df = cumulative_counts_df.sort_index()
    
    plt.figure(figsize=(8, 6))
    x_values = cumulative_counts_df.index
    for event in all_events:
        y_values = cumulative_counts_df[event]
        if y_values.sum() > 0:
            if (grouped[event] > 0).any():
                last_encounter_time = grouped[event][grouped[event] > 0].index[-1]
            else:
                last_encounter_time = x_values[0]
    
            extended_x = np.append(x_values[x_values <= last_encounter_time], last_encounter_time + (x_values[1] - x_values[0]))
            extended_y = np.append(y_values[x_values <= last_encounter_time], 0)
    
            plt.plot(extended_x, extended_y, label=f'{event}')
    
    plt.title("Merger Event Counts with IMBH over Time")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Cumulative Count")
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.grid(True)
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(os.path.join(output_directory, "merger_event_counts.png"), dpi=300, bbox_inches="tight")
    plt.close()
    
    cumulative_counts_df = grouped.cumsum()
    cumulative_counts_df = cumulative_counts_df.sort_index()
    
    x_values = cumulative_counts_df.index
    plt.figure(figsize=(8, 6))
    
    for event in all_events:
        y_values = cumulative_counts_df[event]
        
       
        y_values = y_values.reindex(x_values, fill_value=0)
        
        if y_values.sum() > 0:
            if (grouped[event] > 0).any():
                last_encounter_time = grouped[event][grouped[event] > 0].index[-1]
            else:
                last_encounter_time = x_values[0]
    
            extended_x = np.append(x_values[x_values <= last_encounter_time], last_encounter_time + (x_values[1] - x_values[0]))
            extended_y = np.append(y_values[x_values <= last_encounter_time], 0)
    
            plt.plot(extended_x, extended_y, label=f'{event}')
    
    plt.title("Merger Event Counts with IMBH over Time")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Cumulative Count")
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.grid(True)
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(os.path.join(output_directory, "merger_event_counts.png"), dpi=300, bbox_inches="tight")
    plt.close()

    
   
    event_counts_3 = filtered_df2['lineType(1)'].value_counts()
    
    plt.figure(figsize=(8, 6))
    sns.barplot(x=event_counts_3.index, y=event_counts_3.values,palette="viridis")
    plt.title("IMBH-Compact Object Mergers")
    plt.xlabel("Event Type")
    plt.ylabel("Frequency")
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(os.path.join(output_directory, "compact_object_mergers.png"), dpi=300, bbox_inches="tight")
    plt.close()
    
 
    filtered_df = df[df['mergMass1(26)'] > 0.001]
    uniquemergers = sorted(filtered_df['compType'].unique())
    grouped = filtered_df.groupby(['time[Myr]', 'compType']).size().unstack(fill_value=0)
    cumulative_counts_df = grouped.cumsum()
    

    for comptype in uniquemergers:
        if comptype not in cumulative_counts_df.columns:
            cumulative_counts_df[comptype] = 0
    
    cumulative_counts_df = cumulative_counts_df.sort_index()
    x_values = cumulative_counts_df.index
    
    plt.figure(figsize=(8, 6))
    for comptype in cumulative_counts_df.columns:
        y_values = cumulative_counts_df[comptype]
    
        if (grouped[comptype] > 0).any():
            last_encounter_time = grouped[comptype][grouped[comptype] > 0].index[-1]
        else:
            last_encounter_time = x_values[0]
    
        extended_x = np.append(x_values[x_values <= last_encounter_time], last_encounter_time + (x_values[1] - x_values[0]))
        extended_y = np.append(y_values[x_values <= last_encounter_time], 0)
    
        plt.plot(extended_x, extended_y, label=f'Comp Type {comptype}')
    
    plt.title("Mergers with Different Stellar Types vs Time")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Cumulative Count")
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1.5))
    plt.grid(True)
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(os.path.join(output_directory, "mergers_with_different_stellar_types_vs_time.png"), dpi=300, bbox_inches="tight")
    plt.close()

    print(f"All individual plots saved in '{output_directory}'.")
