import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

file_path = input("Please enter the path to the history file: ")

with open(file_path, 'r') as file:
    history_content = file.readlines()

column_line = history_content[187].strip()

columns = column_line.lstrip('#').split()

if columns[0] == '#':
    columns[0] = 'lineType(1)'

df = pd.read_csv(file_path, delim_whitespace=True, skiprows=188, comment='#', names=columns)

print("DataFrame Head:")
print(df)

print("Column Names:")
print(df.columns.tolist())

#df = df[(df['starType'] == 14)]

df['MassChange'] = df['massNew[Msun]'] - df['massOld[Msun](10)']
mass_increase_points = df[df['MassChange'] > 0.]

sns.set(style="darkgrid")


fig, axs = plt.subplots(3, 2, figsize=(14, 10))
axs.flatten()


axs[0, 0].set_facecolor('black')


axs[0, 0].grid(True, color='white', linestyle='-', linewidth=0.5) 


legend = axs[0, 0].legend(frameon=True)
legend.get_frame().set_facecolor('gray')  
legend.get_frame().set_edgecolor('white')  
plt.setp(legend.get_texts(), color='white')  
axs[0, 0].plot(df['time[Myr]'], df['massNew[Msun]'], label='Mass Evolution', color='w')

ordered_event_types = ['BIN_STAR', 'SIN_EVOL', 'COLLISION', 'BIN_FORM', 'BIN_EVOL', 'BIN_BIN']

marker_styles = {
    'BIN_STAR': 'o',
    'SIN_EVOL': 's',
    'COLLISION': '^',
    'BIN_FORM': 'P',
    'BIN_EVOL': 'X',
    'BIN_BIN': 'D'
}

unique_event_types = df['lineType(1)'].unique()

additional_event_types = [event for event in unique_event_types if event not in ordered_event_types]

remaining_markers = ['v', '<', '>', '*', 'H', '+']  
for i, event in enumerate(additional_event_types):
    marker_styles[event] = remaining_markers[i % len(remaining_markers)]


for event_type in ordered_event_types + additional_event_types:
    event_points = mass_increase_points[mass_increase_points['lineType(1)'] == event_type]
    axs[0, 0].scatter(event_points['time[Myr]'], event_points['massNew[Msun]'],
                      marker=marker_styles[event_type], label=f'{event_type}', s=5, zorder=5)


comp_type_14_points = df[(df['starType'] == 14) & (df['compType'] == 14)]


axs[0, 0].scatter(comp_type_14_points['time[Myr]'], comp_type_14_points['massNew[Msun]'],
                  marker='x', color='yellow', s=10, label='BH-BH', zorder=10)  # s=50 for larger cross marks

axs[0, 0].set_title("Mass Evolution Over Time with Significant Events")
axs[0, 0].set_xlabel("Time (Myr)")
axs[0, 0].set_ylabel("Mass (Msun)")
axs[0, 0].set_xlim([0, 2000])
axs[0, 0].legend(loc='upper left', bbox_to_anchor=(1, 1))

def on_click(event):
    if event.inaxes == axs[0, 0]:  
        # Get the index of the closest point
        closest_index = ((mass_increase_points['time[Myr]'] - event.xdata)**2 + 
                         (mass_increase_points['massNew[Msun]'] - event.ydata)**2).idxmin()

        
        event_info = mass_increase_points.loc[closest_index]

        output_text = f"""
        Event Information:
        --------------------
        EventType: {event_info['lineType(1)']}
        EncId: {event_info['encId']}
        EncCase: {event_info['encCase']}
        Time (Myr): {event_info['time[Myr]']}
        Summary: {event_info['Summary']}
        Id6: {event_info['id(6)']}
        IdComp: {event_info['idComp']}
        IdCompNew: {event_info['idCompNew']}
        IM: {event_info['im']}
        MassOld (Msun): {event_info['massOld[Msun](10)']}
        MassNew (Msun): {event_info['massNew[Msun]']}
        MassChange: {event_info['MassChange']}
        aOld (Rsun): {event_info['aOld[Rsun]']}
        aNew (Rsun): {event_info['aNew[Rsun]']}
        eOld: {event_info['eOld(16)']}
        eNew: {event_info['eNew']}
        starType: {event_info['starType']}
        compType: {event_info['compType']}
        binType: {event_info['binType']}
        
        """
        print(output_text)



fig.canvas.mpl_connect('button_press_event', on_click)



event_counts = df['lineType(1)'].value_counts()
sns.barplot(x=event_counts.index, y=event_counts.values, ax=axs[0, 1])
axs[0, 1].set_title("Event Type Frequency")
axs[0, 1].set_xlabel("Event Type")
axs[0, 1].set_ylabel("Frequency")
axs[0, 1].tick_params(axis='x', rotation=45)


axs[1, 0].plot(df['time[Myr]'], df['r[pc]'], label='r', color='r')
axs[1, 0].set_title("Position Evolution Over Time")
axs[1, 0].set_xlabel("Time (Myr)")
axs[1, 0].set_ylabel("r [pc]")
axs[1, 0].legend()


axs[1, 1].plot(df['time[Myr]'], df['eNew'], label='eNew', color='r')
axs[1, 1].set_title("Orbital Eccentricity Evolution Over Time")
axs[1, 1].set_xlabel("Time (Myr)")
axs[1, 1].set_ylabel("Orbital Eccentricity")
axs[1, 1].legend()


filtered_df = df[df['encCase'] != 0]

unique_enc_cases = sorted(filtered_df['encCase'].unique())

grouped = filtered_df.groupby(['time[Myr]', 'encCase']).size().unstack(fill_value=0)


cumulative_counts_df = grouped.cumsum()


for enc_case in unique_enc_cases:
    if enc_case not in cumulative_counts_df.columns:
        cumulative_counts_df[enc_case] = 0


cumulative_counts_df = cumulative_counts_df.sort_index()

x_values = cumulative_counts_df.index  

for enc_case in cumulative_counts_df.columns:
    y_values = cumulative_counts_df[enc_case]
    
  
    if (grouped[enc_case] > 0).any():
        last_encounter_time = grouped[enc_case][grouped[enc_case] > 0].index[-1]  
    else:
        last_encounter_time = x_values[0]  

    
    extended_x = np.append(x_values[x_values <= last_encounter_time], last_encounter_time + (x_values[1] - x_values[0]))  # Add a point after the last occurrence
    extended_y = np.append(y_values[x_values <= last_encounter_time], 0) 
    
    
    axs[2, 0].plot(extended_x, extended_y, label=f'Encounter Case {enc_case}')


axs[2, 0].set_title("Cumulative Encounter Cases Over Time")
axs[2, 0].set_xlabel("Time (Myr)")
axs[2, 0].set_ylabel("Cumulative Count")
axs[2, 0].legend(loc='upper left', bbox_to_anchor=(1, 1))
axs[2, 0].grid(True)
axs[2,0].set_yscale('log')
plt.tight_layout()




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


blackhole_interactions_df = filtered_df[(filtered_df['starType'] == 14) & (filtered_df['compType'] == 14)]

mergers_with_blackholes = 0
flybys_with_blackholes = 0
exchanges_with_blackholes = 0
single_stars_with_blackholes = 0

for enc_case in blackhole_interactions_df['encCase']:
    enc_case_str = str(enc_case).zfill(3)  
   
    first_digit = int(enc_case_str[0])
    if first_digit in [1, 3, 6]:
        mergers_with_blackholes += 1  
    
    # Second digit: Count flybys/exchanges
    second_digit = int(enc_case_str[1])
    if second_digit in [1, 2]:
        flybys_with_blackholes += 1  
    elif second_digit == 3:
        exchanges_with_blackholes += 1  #single exchange
    elif second_digit == 6:
        exchanges_with_blackholes += 1  #double exchange
    
  
    third_digit = int(enc_case_str[2])
    single_stars_with_blackholes += third_digit  

axs[2, 1].axis('off')  
stats_text = f"""

Mergers:
- Same binary mergers: {mergers_same_binary}
- Different binary mergers: {mergers_diff_binary}
- Double mergers with exchanges: {double_mergers}
- Mergers with blackholes: {mergers_with_blackholes}

Flybys/Exchanges:
- No change in binaries: {binaries_no_change}
- Single exchange: {exchanges_single}
- Double exchange: {exchanges_double}
- Binary unbounded or merged: {no_binaries_left}

Flybys/Exchanges with Blackholes:
- Flybys with blackholes: {flybys_with_blackholes}
- Single exchange with blackholes: {exchanges_with_blackholes}
- Double exchange with blackholes: {exchanges_with_blackholes}

Single Stars formed after interactions:
- Total: {single_stars_count}
- Single stars formed with blackholes: {single_stars_with_blackholes}
"""


axs[2, 1].text(0.5, 0.5, stats_text, ha='center', va='center', fontsize=10, bbox=dict(facecolor='lightgray', edgecolor='black'))
#axs[2, 1].set_title("Encounter Stats")



plt.tight_layout()
plt.show()


fig, axs = plt.subplots(3, 2, figsize=(14, 10))
axs.flatten()

filtered_df1 = df[df['mergMass1(26)'] > 0.001]

event_counts_2 = filtered_df1['lineType(1)'].value_counts()
sns.barplot(x=event_counts_2.index, y=event_counts_2.values, ax=axs[0, 0])
axs[0, 0].set_title("Mergers")
axs[0, 0].set_xlabel("Event Type")
axs[0, 0].set_ylabel("Frequency")
axs[0, 0].tick_params(axis='x', rotation=45)

filtered_df2 = df[df['mergMass1(26)'] > 0.001]

event_counts_3 = filtered_df2['compType'].value_counts()
sns.barplot(x=event_counts_3.index, y=event_counts_3.values, ax=axs[0, 1])
axs[0, 1].set_title("Mergers with each stellar type")
axs[0, 1].set_xlabel("Stellar Type")
axs[0, 1].set_ylabel("Frequency")


filtered_df3 = filtered_df2[(filtered_df2['starType'] ==14) & (filtered_df2['compType']==14)]
filtered_df_gw = filtered_df3[filtered_df3['lineType(1)'] == 'BIN_EVOL']

bh1= df[df['idComp']==204]
axs[1, 0].plot(filtered_df_gw['time[Myr]'], filtered_df_gw['aNew[Rsun]'])
axs[1, 0].set_title("Evolution of semi-major axis in BH-BH Binary Mergers")

filtered_df4=df[df['lineType(1)'] == 'BIN_FORM']

grouped_bin_form = filtered_df4.groupby('time[Myr]').size()


cumulative_bin_form = grouped_bin_form.cumsum()

axs[1, 1].plot(cumulative_bin_form.index, cumulative_bin_form.values, label='Cumulative BIN_FORM Events', color='blue', linewidth=2)
axs[1, 1].set_title("New binary formation over time")


filtered_df = df[df['mergMass1(26)'] > 0.001]
lowmass=filtered_df[filtered_df['mergStarType1']==0]
print(lowmass.head())
print('low mass main sequence in mergers: ', len(lowmass))

uniquemergers = sorted(filtered_df['mergStarType1'].unique())


grouped = filtered_df.groupby(['time[Myr]', 'mergStarType1']).size().unstack(fill_value=0)


cumulative_counts_df = grouped.cumsum()


for comptype in uniquemergers:
    if comptype not in cumulative_counts_df.columns:
        cumulative_counts_df[comptype] = 0


cumulative_counts_df = cumulative_counts_df.sort_index()


x_values = cumulative_counts_df.index


for comptype in cumulative_counts_df.columns:
    y_values = cumulative_counts_df[comptype]
    
   
    if (grouped[comptype] > 0).any():
        last_encounter_time = grouped[comptype][grouped[comptype] > 0].index[-1]
    else:
        last_encounter_time = x_values[0]  


    extended_x = np.append(x_values[x_values <= last_encounter_time], last_encounter_time + (x_values[1] - x_values[0]))
    extended_y = np.append(y_values[x_values <= last_encounter_time], 0)  

  
    axs[2, 0].plot(extended_x, extended_y, label=f'Comp Type {comptype}')


axs[2, 0].set_title("Encounters with Different Stellar Types")
axs[2, 0].set_xlabel("Time (Myr)")
axs[2, 0].set_ylabel("Cumulative Count")
axs[2, 0].legend(loc='upper left', bbox_to_anchor=(1, 1))
axs[2, 0].grid(True)
axs[2, 0].set_yscale('log')

plt.tight_layout()
plt.show()