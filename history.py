import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import random
import os
import warnings
import config
import logging
import os



file_path = input("Please enter the path to the history file: ")

history_basename = os.path.splitext(os.path.basename(file_path))[0]

with open(file_path, 'r') as file:
    history_content = file.readlines()
    
    column_line = history_content[187].strip()
    
    columns = column_line.lstrip('#').split()
    
    if columns[0] == '#':
        columns[0] = 'lineType(1)'

df = pd.read_csv(file_path, sep='\s+', skiprows=188, comment='#', names=columns)

print("DataFrame Head:")
print(df)

print("Column Names:")
print(df.columns.tolist())

def save_data(df, filename):

    output_directory = f"data_{history_basename}"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    file_path = os.path.join(output_directory, filename)
    
    
    column_widths = []
    for col in df.columns:
        max_col_width = max(df[col].astype(str).apply(len).max(), len(col))
        column_widths.append(max_col_width)
    
    def format_row(row):
        return '\t\t'.join(f'{str(item):<{width}}' for item, width in zip(row, column_widths))
    with open(file_path, 'w') as f:
        f.write(format_row(df.columns) + '\n')
        
        for index, row in df.iterrows():
            f.write(format_row(row.values) + '\n')

print(f"Data saved to {file_path}")
df = df[(df['starType'] == 14)]

df['MassChange'] = df['massNew[Msun]'] - df['massOld[Msun](10)']
mass_increase_points = df[df['MassChange'] > 0.]




# Filter and calculate mass change
df = df[df['starType'] == 14]
df['MassChange'] = df['massNew[Msun]'] - df['massOld[Msun](10)']
mass_increase_points = df[df['MassChange'] > 0.]

import matplotlib.pyplot as plt
import seaborn as sns

# Set publication-friendly style
sns.set(style='whitegrid', context='paper', font_scale=1.2)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['axes.linewidth'] = 1.2

fig, ax = plt.subplots(figsize=(10, 6))

# Plot the mass evolution curve and significant event markers
ax.plot(df['time[Myr]'], df['massNew[Msun]'], color='red', linewidth=2, label='Mass Evolution')
ax.scatter(mass_increase_points['time[Myr]'], mass_increase_points['massNew[Msun]'],
           color='blue', marker='o', s=50, zorder=3, label='Significant Event')

# Titles, labels, and legend
ax.set_title("Mass Evolution Over Time with Significant Events", fontsize=14, fontweight='bold')
ax.set_xlabel("Time (Myr)", fontsize=12)
ax.set_ylabel("Mass (Msun)", fontsize=12)
ax.legend(loc='upper right', fontsize=10)

# Let the axes autoscale with a small margin
ax.autoscale(enable=True, axis='both', tight=True)
ax.margins(0.1)
ax.grid(True, linestyle='--', linewidth=0.7, color='gray', alpha=0.7)

plt.tight_layout()
plt.savefig('mass_evolution_improved.png', dpi=300)
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style='whitegrid', context='paper', font_scale=1.2)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['axes.linewidth'] = 1.2

fig, ax = plt.subplots(figsize=(10, 6))

# Full mass evolution curve
ax.plot(df['time[Myr]'], df['massNew[Msun]'],
        color='red', linewidth=2, label='Mass Evolution')

# Scatter plot for all significant events (no time limit)
ax.scatter(mass_increase_points['time[Myr]'], mass_increase_points['massNew[Msun]'],
           color='blue', marker='o', s=50, zorder=3, label='Significant Event')

ax.set_title("Mass Evolution Over Time with Significant Events (Entire Evolution)",
             fontsize=14, fontweight='bold')
ax.set_xlabel("Time (Myr)", fontsize=12)
ax.set_ylabel("Mass (Msun)", fontsize=12)
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1))
ax.grid(True, linestyle='--', linewidth=0.7, color='gray', alpha=0.7)

plt.tight_layout()
plt.savefig('mass_evolution_full.png', dpi=300)
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns

fig, ax = plt.subplots(figsize=(8, 6))

# Compute event counts and define a color palette
event_counts = df['lineType(1)'].value_counts()
palette = sns.color_palette("Set1", n_colors=len(event_counts))

sns.barplot(x=event_counts.index, y=event_counts.values, palette=palette, ax=ax)
ax.set_title("Event Type Frequency", fontsize=14, fontweight='bold')
ax.set_xlabel("Event Type", fontsize=12)
ax.set_ylabel("Frequency", fontsize=12)
ax.tick_params(axis='x', rotation=45)
ax.grid(True, linestyle='--', linewidth=0.7, alpha=0.7)

plt.tight_layout()
plt.savefig('event_type_frequency.png', dpi=300)
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(df['time[Myr]'], df['r[pc]'], color='red', linewidth=2, label='IMBH Position')
ax.set_title("IMBH Position Change Over Time", fontsize=14, fontweight='bold')
ax.set_xlabel("Time (Myr)", fontsize=12)
ax.set_ylabel("r [pc]", fontsize=12)
ax.set_xscale('log')
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1))
ax.grid(True, linestyle='--', linewidth=0.7, alpha=0.7)

plt.tight_layout()
plt.savefig('imbh_position_change.png', dpi=300)
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

fig, ax = plt.subplots(figsize=(10, 6))

# Filter out encounters and group by time and encounter case
filtered_df = df[df['encCase'] != 0]
unique_enc_cases = sorted(filtered_df['encCase'].unique())
grouped = filtered_df.groupby(['time[Myr]', 'encCase']).size().unstack(fill_value=0)
cumulative_counts_df = grouped.cumsum().sort_index()
x_values = cumulative_counts_df.index

# Plot cumulative counts for each encounter case
for enc_case in unique_enc_cases:
    y_values = cumulative_counts_df[enc_case]
    # Determine last time where this encounter case appears
    if (grouped[enc_case] > 0).any():
        last_encounter_time = grouped[enc_case][grouped[enc_case] > 0].index[-1]
    else:
        last_encounter_time = x_values[0]
    step = x_values[1] - x_values[0] if len(x_values) > 1 else 1
    extended_x = np.append(x_values[x_values <= last_encounter_time],
                           last_encounter_time + step)
    extended_y = np.append(y_values[x_values <= last_encounter_time], 0)
    ax.plot(extended_x, extended_y, label=f'Encounter Case {enc_case}', linewidth=2)

ax.set_title("Cumulative Encounter Cases Over Time", fontsize=14, fontweight='bold')
ax.set_xlabel("Time (Myr)", fontsize=12)
ax.set_ylabel("Cumulative Count", fontsize=12)
ax.set_yscale('log')
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1))
ax.grid(True, linestyle='--', linewidth=0.7, alpha=0.7)

plt.tight_layout()
plt.savefig('cumulative_encounter_cases.png', dpi=300)
plt.show()


import matplotlib.pyplot as plt

mergers_same_binary = 0
mergers_diff_binary = 0
double_mergers = 0
binaries_no_change = 0
exchanges_single = 0
exchanges_double = 0
no_binaries_left = 0
single_stars_count = 0
 
 
for enc_case in filtered_df['encCase']:
    enc_case_str = str(enc_case).zfill(3)  
    # First digit analysis (mergers)
    first_digit = int(enc_case_str[0])
    if first_digit == 1:
        mergers_same_binary += 1
    elif first_digit == 3:
        mergers_diff_binary += 1
    elif first_digit == 6:
        double_mergers += 1

    # Second digit analysis (flybys/exchanges)
    second_digit = int(enc_case_str[1])
    if second_digit in [1, 2]:
        binaries_no_change += 1
    elif second_digit == 3:
        exchanges_single += 1
    elif second_digit == 6:
        exchanges_double += 1
    elif second_digit == 0:
        no_binaries_left += 1
 
    # Third digit analysis (single stars)
    third_digit = int(enc_case_str[2])
    single_stars_count += third_digit
fig, ax = plt.subplots(figsize=(8, 6))
ax.axis('off')

# Assume these variables have been computed earlier:
# mergers_same_binary, mergers_diff_binary, double_mergers, binaries_no_change,
# exchanges_single, exchanges_double, no_binaries_left, single_stars_count

stats_text = f"""
Mergers:
- Same binary mergers: {mergers_same_binary}
- Different binary mergers: {mergers_diff_binary}
- Double mergers with exchanges: {double_mergers}

Flybys/Exchanges:
- No change in binaries: {binaries_no_change}
- Single exchange: {exchanges_single}
- Double exchange: {exchanges_double}
- Binary unbounded/merged: {no_binaries_left}

Single Stars formed after interactions:
- Total: {single_stars_count}
"""

ax.text(0.5, 0.5, stats_text, ha='center', va='center', fontsize=12,
        bbox=dict(facecolor='lightgray', edgecolor='black', pad=10))
plt.tight_layout()
plt.savefig('encounter_stats.png', dpi=300)
plt.show()

import matplotlib.pyplot as plt
import seaborn as sns

fig, ax = plt.subplots(figsize=(8, 6))

# Filter mergers events (assume mergMass1(26) > 0.001 indicates a merger)
filtered_df1 = df[df['mergMass1(26)'] > 0.001]
event_counts_2 = filtered_df1['lineType(1)'].value_counts()

sns.barplot(x=event_counts_2.index, y=event_counts_2.values, palette="Set1", ax=ax)
ax.set_title("Mergers Events", fontsize=14, fontweight='bold')
ax.set_xlabel("Event Type", fontsize=12)
ax.set_ylabel("Frequency (log scale)", fontsize=12)
ax.set_yscale('log')
ax.tick_params(axis='x', rotation=45)
ax.grid(True, linestyle='--', linewidth=0.7, alpha=0.7)

plt.tight_layout()
plt.savefig('mergers_events.png', dpi=300)
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns

fig, ax = plt.subplots(figsize=(8, 6))

filtered_df2 = df[df['mergMass1(26)'] > 0.001].copy()
event_counts_3 = filtered_df2['compType'].value_counts()

sns.barplot(x=event_counts_3.index, y=event_counts_3.values, palette="Set1", ax=ax)
ax.set_title("Mergers with Each Stellar Type", fontsize=14, fontweight='bold')
ax.set_xlabel("Stellar Type", fontsize=12)
ax.set_ylabel("Frequency (log scale)", fontsize=12)
ax.set_yscale('log')
ax.tick_params(axis='x', rotation=45)
ax.grid(True, linestyle='--', linewidth=0.7, alpha=0.7)

plt.tight_layout()
plt.savefig('mergers_stellar_type.png', dpi=300)
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

fig, ax = plt.subplots(figsize=(10, 6))

filtered_df2 = df[df['mergMass1(26)'] > 0.001].copy()
all_events = sorted(filtered_df2['lineType(1)'].unique())
grouped = filtered_df2.groupby(['time[Myr]', 'lineType(1)']).size().unstack(fill_value=0)
cumulative_counts_df = grouped.cumsum().sort_index()
x_values = cumulative_counts_df.index

for event in all_events:
    y_values = cumulative_counts_df[event]
    if (grouped[event] > 0).any():
        last_encounter_time = grouped[event][grouped[event] > 0].index[-1]
    else:
        last_encounter_time = x_values[0]
    step = x_values[1] - x_values[0] if len(x_values) > 1 else 1
    extended_x = np.append(x_values[x_values <= last_encounter_time],
                           last_encounter_time + step)
    extended_y = np.append(y_values[x_values <= last_encounter_time], 0)
    ax.plot(extended_x, extended_y, label=f'{event}', linewidth=2)

ax.set_title("Merger Event Counts with IMBH Over Time", fontsize=14, fontweight='bold')
ax.set_xlabel("Time (Myr)", fontsize=12)
ax.set_ylabel("Cumulative Count (log scale)", fontsize=12)
ax.set_yscale('log')
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1))
ax.grid(True, linestyle='--', linewidth=0.7, alpha=0.7)

plt.tight_layout()
plt.savefig('merger_event_counts.png', dpi=300)
plt.show()


import matplotlib.pyplot as plt
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
        count = len(filtered_df_by_type[filtered_df_by_type['lineType(1)'] == line_type])
        summary_table[line_type].append(count)
        total_count += count
    
    summary_table['Total Count'].append(total_count)
    
summary_df = pd.DataFrame(summary_table)

fig, ax = plt.subplots(figsize=(8, 6))
ax.axis('off')

# Assume summary_df is your summary DataFrame computed from filtered_df2
table = ax.table(cellText=summary_df.values, colLabels=summary_df.columns,
                 cellLoc='center', loc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.5, 1.5)

# Format header and alternate row colors for readability
for (row, col), cell in table.get_celld().items():
    if row == 0:
        cell.set_text_props(weight='bold', color='white')
        cell.set_facecolor('#40466e')
    else:
        cell.set_facecolor('#f1f1f2' if row % 2 == 0 else 'white')

plt.tight_layout()
plt.savefig('summary_table.png', dpi=300)
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

fig, ax = plt.subplots(figsize=(10, 6))

filtered_df = filtered_df2.copy()
uniquemergers = sorted(filtered_df['compType'].unique())
grouped = filtered_df.groupby(['time[Myr]', 'compType']).size().unstack(fill_value=0)
cumulative_counts_df = grouped.cumsum().sort_index()
x_values = cumulative_counts_df.index

# Define a color palette for consistency
color_palette = [
    'crimson', 'royalblue', 'goldenrod', 'mediumseagreen',
    'tomato', 'slateblue', 'darkorange', 'mediumorchid',
    'lightcoral', 'seagreen', 'slategray', 'peru', 'steelblue', 'sandybrown'
]

for i, comptype in enumerate(cumulative_counts_df.columns):
    y_values = cumulative_counts_df[comptype]
    if (grouped[comptype] > 0).any():
        last_encounter_time = grouped[comptype][grouped[comptype] > 0].index[-1]
    else:
        last_encounter_time = x_values[0]
    step = (x_values[1] - x_values[0]) if len(x_values) > 1 else 1
    extended_x = np.append(x_values[x_values <= last_encounter_time],
                           last_encounter_time + step)
    extended_y = np.append(y_values[x_values <= last_encounter_time], 0)
    ax.plot(extended_x, extended_y, label=f'Comp Type {comptype}', linewidth=2,
            color=color_palette[i % len(color_palette)])

ax.set_title("Mergers with Different Stellar Types vs. Time", fontsize=14, fontweight='bold')
ax.set_xlabel("Time (Myr)", fontsize=12)
ax.set_ylabel("Cumulative Count (log scale)", fontsize=12)
ax.set_yscale('log')
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.5))
ax.grid(True, linestyle='--', linewidth=0.7, alpha=0.7)

plt.tight_layout()
plt.savefig('mergers_stellar_types_vs_time.png', dpi=300)
plt.show()

import matplotlib.pyplot as plt
import seaborn as sns

fig, ax = plt.subplots(figsize=(8, 6))

filtered_df3 = filtered_df2[(filtered_df2['starType'] == 14) & (filtered_df2['compType'].isin([10,11,12,13, 14]))]
filtered_df_gw = filtered_df3[filtered_df3['lineType(1)'] == 'BIN_EVOL']
# Assume filtered_df3 is defined as a subset of filtered_df2 with specific conditions
event_counts_4 = filtered_df3['lineType(1)'].value_counts()
palette = sns.color_palette("Set1", n_colors=len(event_counts_4.index))

sns.barplot(x=event_counts_4.index, y=event_counts_4.values, palette=palette, ax=ax)
ax.set_title("IMBH–Compact Object Mergers", fontsize=14, fontweight='bold')
ax.set_xlabel("Event Type", fontsize=12)
ax.set_ylabel("Frequency (log scale)", fontsize=12)
ax.set_yscale('log')
ax.tick_params(axis='x', rotation=45)
ax.grid(True, linestyle='--', linewidth=0.7, alpha=0.7)

plt.tight_layout()
plt.savefig('imbh_compact_object_mergers.png', dpi=300)
plt.show()




