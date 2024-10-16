#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 14:55:02 2024

@author: sohaib
"""
import glob
import pandas as pd
from sanitization import validate_file_path
from utils import extract_file_name, save_data


try:
    while True:
        file_path = input("Enter path to history files (with extension): ")
        try:
            validate_file_path(file_path)
            break
        except ValueError as e:
            print(e)
except Exception as e:
    print(e)


for file_path in glob.glob(file_path):
    file_name = extract_file_name(file_path)
    print(f"Processing file: {file_name}")
    history_basename = file_name.split('.')[0]

    with open(file_path, 'r') as file:
        history_content = file.readlines()

    column_line = history_content[86].strip()

    columns = column_line.lstrip('#').split()

    if columns[0] == '#':
        columns[0] = 'lineType(1)'

    df = pd.read_csv(file_path, delim_whitespace=True,
                     skiprows=87, comment='#', names=columns)

    print("DataFrame Head:")
    print(df)

    print("Column Names:")
    print(df.columns.tolist())

    df = df[(df['starType'] == 14)]

    df['MassChange'] = df['massNew[Msun](10)'] - df['massOld[Msun]']
    mass_increase_points = df[df['MassChange'] > 0.]

    comp_type_14_points = df[(df['starType'] == 14) & (df['compType'] == 14)]

    filtered_df1 = df[df['mergMass1(25)'] > 0.001]

    filtered_df2 = df[df['mergMass1(25)'] > 0.001]
    save_data(filtered_df2, 'mergers.dat', history_basename)

    filtered_df3 = filtered_df2[(filtered_df2['starType'] == 14) & (
        filtered_df2['compType'].isin([10, 11, 12, 13, 14]))]
    filtered_df_gw = filtered_df3[filtered_df3['lineType(1)'] == 'BIN_EVOL']
    print("Gravitational waves candidates: ", filtered_df_gw)
    save_data(filtered_df_gw, 'gwcandidates.dat', history_basename)

    if not filtered_df_gw.empty:

        unique_ids_gw = filtered_df_gw['idComp'].unique()

        for id in unique_ids_gw:
            gw = df[df['idComp'] == id]

            if not gw.empty:
                comp_type = gw['compType'].iloc[0]
                if comp_type in [10, 11, 12]:
                    filename = f'gwcandidate_bh-wd{id}.dat'
                    save_data(gw, filename, history_basename)
                elif comp_type == 13:
                    filename = f'gwcandidate_bh-ns{id}.dat'
                    save_data(gw, filename, history_basename)
                elif comp_type == 14:
                    filename = f'gwcandidate_bh-bh{id}.dat'
                    save_data(gw, filename, history_basename)

    print("History files for GW candidates have been generated")

    unique_line_types = filtered_df2['lineType(1)'].unique()

    summary_table = {
        'Stellar Type': [],
    }

    for line_type in unique_line_types:
        summary_table[line_type] = []

    summary_table['Total Count'] = []

    for comp_type in range(15):
        filtered_df_by_type = filtered_df2[filtered_df2['compType'] == comp_type]

        summary_table['Stellar Type'].append(comp_type)

        total_count = 0
        for line_type in unique_line_types:
            count = len(
                filtered_df_by_type[filtered_df_by_type['lineType(1)'] == line_type])
            summary_table[line_type].append(count)
            total_count += count

        summary_table['Total Count'].append(total_count)

    summary_df = pd.DataFrame(summary_table)

    save_data(summary_df, 'mergersum.dat', history_basename)

    del df

print("Finished processing all files.")
