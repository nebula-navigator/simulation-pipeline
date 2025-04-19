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



logging.getLogger().setLevel(logging.CRITICAL)




def history():
    
    

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
    #df = df[(df['starType'] == 14)]
    
    

    
    
        # Filter and calculate mass change
   # df = df[df['starType'] == 14]
    df['MassChange'] = df['massNew[Msun]'] - df['massOld[Msun](10)']
    mass_increase_points = df[df['MassChange'] != 0.]
    
 
    
  
    fig, axs = plt.subplots(3, 2, figsize=(14, 10))
    
   
    ax_zoom = axs[0, 0]
    ax_full = axs[1, 1]
    

    for ax in (ax_zoom, ax_full):
        ax.set_facecolor('white')
        ax.grid(True, color='black', linestyle='-', linewidth=0.5)
    
    
    ax_zoom.plot(df['time[Myr]'], df['massNew[Msun]'],
                 label='Mass Evolution', color='r')
    ax_full.plot(df['time[Myr]'], df['massNew[Msun]'],
                 label='Mass Evolution', color='r')
    
   
        
    legend1 = ax_zoom.legend(frameon=True)
    legend1.get_frame().set_facecolor('white')  
    legend1.get_frame().set_edgecolor('gray')
    plt.setp(legend1.get_texts(), color='black')
    
    legend2 = ax_full.legend(frameon=True)
    legend2.get_frame().set_facecolor('white')
    legend2.get_frame().set_edgecolor('gray')
    plt.setp(legend2.get_texts(), color='black')
    
  
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
        ax_zoom.scatter(event_points['time[Myr]'], event_points['massNew[Msun]'],
                        marker=marker_styles[event_type], label=f'{event_type}', s=4, zorder=2)
        ax_full.scatter(event_points['time[Myr]'], event_points['massNew[Msun]'],
                        marker=marker_styles[event_type], label=f'{event_type}', s=4, zorder=2)
    
   
    # comp_type_14_points = df[(df['starType'] == 14) & (df['compType'] == 14)]
    # ax_zoom.scatter(comp_type_14_points['time[Myr]'], comp_type_14_points['massNew[Msun]'],
    #                 marker='x', color='yellow', s=10, label='BH-BH', zorder=10)
    # ax_full.scatter(comp_type_14_points['time[Myr]'], comp_type_14_points['massNew[Msun]'],
    #                 marker='x', color='yellow', s=10, label='BH-BH', zorder=10)
    
   
    ax_zoom.set_title("Mass Evolution Over Time with Significant Events (0-2Gyr)")
    ax_zoom.set_xlabel("Time (Myr)")
    ax_zoom.set_ylabel("Mass (Msun)")
    ax_zoom.set_xlim([0, 2000])
    ax_zoom.legend(loc='upper left', bbox_to_anchor=(1, 1))
    
    ax_full.set_title("Mass Evolution Over Time with Significant Events (Entire Evolution)")
    ax_full.set_xlabel("Time (Myr)")
    ax_full.set_ylabel("Mass (Msun)")
    ax_full.legend(loc='upper left', bbox_to_anchor=(1, 1))
    

    def on_click(event):
        if event.inaxes == ax_zoom or event.inaxes == ax_full:
            # Find the index of the closest point
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
    
    
    import matplotlib.patches as mpatches
   

    event_counts = df['lineType(1)'].value_counts()
    palette = sns.color_palette("Set1", n_colors=len(event_counts.index))
    sns.barplot(
        x=event_counts.index,
        y=event_counts.values,
        ax=axs[0, 1]
        # hue=event_counts.index,   
        # palette=palette,
        # dodge=False,
        # legend=False
    )
    axs[0, 1].set_title("Event Type Frequency")
    axs[0, 1].set_xlabel("Event Type")
    axs[0, 1].set_ylabel("Frequency")
    axs[0, 1].set_yscale('log')
    axs[0, 1].tick_params(axis='x', rotation=45)
    handles = [
        mpatches.Patch(color=palette[i], label=str(count))
        for i, count in enumerate(event_counts.values)
    ]
    #axs[0, 1].legend(handles=handles, loc='upper left', bbox_to_anchor=(1, 1))

  
    
    axs[1, 0].plot(df['time[Myr]'], df['r[pc]'], label='r', color='r')
    axs[1, 0].set_title("IMBH Position Change Over Time")
    axs[1, 0].set_xlabel("Time (Myr)")
    axs[1, 0].set_xscale('log')
    axs[1, 0].set_ylabel("r [pc]")
    axs[1, 0].legend()
    
    
    # axs[1, 1].plot(df['time[Myr]'], df['eNew'], label='eNew', color='r')
    # axs[1, 1].set_title("Orbital Eccentricity Evolution Over Time")
    # axs[1, 1].set_xlabel("Time (Myr)")
    # axs[1, 1].set_ylabel("Orbital Eccentricity")
    # axs[1, 1].legend()
    
    
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
    
        
        extended_x = np.append(x_values[x_values <= last_encounter_time], last_encounter_time + (x_values[1] - x_values[0]))  
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
    
    
     
    
    axs[2, 1].axis('off')  
    stats_text = f"""
    
    Mergers:
    - Same binary mergers: {mergers_same_binary}
    - Different binary mergers: {mergers_diff_binary}
    - Double mergers with exchanges: {double_mergers}
    
    
    Flybys/Exchanges:
    - No change in binaries: {binaries_no_change}
    - Single exchange: {exchanges_single}
    - Double exchange: {exchanges_double}
    - Binary unbounded or merged: {no_binaries_left}
    
    
    
    Single Stars formed after interactions:
    - Total: {single_stars_count}
    
    """
    
    
    axs[2, 1].text(0.5, 0.5, stats_text, ha='center', va='center', fontsize=10, bbox=dict(facecolor='lightgray', edgecolor='black'))
    #axs[2, 1].set_title("Encounter Stats")
    
    
    plt.suptitle("IMBH MASS EVOLUTON ANALYSIS", fontsize=16)
    
    plt.tight_layout()
    
    plt.savefig('subplots1.png', dpi=300)
    

    
    fig, axs = plt.subplots(3, 2, figsize=(14, 10))
    axs.flatten()
    
    filtered_df1 = df[df['mergMass1(26)'] > 0.001]
    
    event_counts_2 = filtered_df1['lineType(1)'].value_counts()
    x=event_counts_2
    hue=x
    #sns.barplot(x=event_counts_2.index, y=event_counts_2.values, ax=axs[0, 0],palette="Set1",hue=hue)
    sns.barplot(x=event_counts_2.index, y=event_counts_2.values, ax=axs[0, 0])
    axs[0, 0].set_title("Mergers Events")
    axs[0, 0].set_xlabel("Event Type")
    axs[0, 0].set_ylabel("Frequency")
    axs[0, 0].set_yscale('log')
    axs[0, 0].tick_params(axis='x')
   # axs[0, 0].legend(loc='upper left', bbox_to_anchor=(1, 1))
    
    
    
    filtered_df2 = df[df['mergMass1(26)'] > 0.001].copy()
    save_data(filtered_df2,'mergers.dat')
    
    event_counts_3 = filtered_df2['compType'].value_counts()
    x=event_counts_3
    hue=x
    color_palette = [
        'crimson',       # #DC143C
        'royalblue',     # #4169E1
        'goldenrod',     # #DAA520
        'mediumseagreen',# #3CB371
        'tomato',        # #FF6347
        'slateblue',     # #6A5ACD
        'darkorange',    # #FF8C00
        'mediumorchid',  # #BA55D3
        'lightcoral',    # #F08080
        'seagreen',      # #2E8B57
        'slategray',     # #708090
        'peru',          # #CD853F
        'steelblue',     # #4682B4
        'sandybrown'     # #F4A460
    ]
    
    #sns.barplot(x=event_counts_3.index, y=event_counts_3.values, ax=axs[0, 1],palette = color_palette,hue=hue)
    sns.barplot(x=event_counts_3.index, y=event_counts_3.values, ax=axs[0, 1])
    axs[0, 1].set_title("Mergers with each stellar type")
    axs[0, 1].set_xlabel("Stellar Type")
    axs[0, 1].set_ylabel("Frequency")
    axs[0, 1].tick_params(axis='x')
    axs[0, 1].set_yscale('log')
    #axs[0, 1].legend(loc='upper left', bbox_to_anchor=(1, 1.2))
    
    filtered_df3 = filtered_df2[(filtered_df2['starType'] == 14) & (filtered_df2['compType'].isin([10,11,12,13, 14]))]
    filtered_df_gw = filtered_df3[filtered_df3['lineType(1)'] == 'BIN_EVOL']
    print("Gravitational waves candidates: ", filtered_df_gw)
    save_data(filtered_df_gw,'gwcandidates.dat')
    
    unique_ids_gw = filtered_df_gw['idComp'].unique()
    
    for id in unique_ids_gw:
        gw = df[df['idComp'] == id]
    
        if not gw.empty:
            comp_type = gw['compType'].iloc[0]  
            if comp_type in [10, 11, 12]:
                filename = f'gwcandidate_bh-wd{id}.dat'
                save_data(gw, filename)
            elif comp_type == 13:
                filename = f'gwcandidate_bh-ns{id}.dat'
                save_data(gw, filename)
            elif comp_type == 14:
                filename = f'gwcandidate_bh-bh{id}.dat'
                save_data(gw, filename)
    
    print("History files for GW candidates have been generated")
    
        
    filtered_df4 = filtered_df2.copy()
    all_events = sorted(filtered_df2['lineType(1)'].unique())
    grouped = filtered_df2.groupby(['time[Myr]', 'lineType(1)']).size().unstack(fill_value=0).reindex(columns=all_events, fill_value=0)
    
    
    cumulative_counts_df = grouped.cumsum()
    cumulative_counts_df = cumulative_counts_df.sort_index()
    
    
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
            
            axs[1, 0].plot(extended_x, extended_y, label=f'{event}')


    
    axs[1, 0].set_title("Merger event counts with IMBH over Time")
    axs[1, 0].set_xlabel("Time (Myr)")
    axs[1, 0].set_ylabel("Cumulative Count")
    axs[1, 0].legend(loc='upper left', bbox_to_anchor=(1, 1))
    axs[1, 0].grid(True)
    axs[1, 0].set_yscale('log')
    plt.tight_layout()
    

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
    
    axs[1, 1].axis('off')
    
    table_data = summary_df.values
    column_labels = summary_df.columns
    
    table = axs[1, 1].table(
        cellText=table_data,
        colLabels=column_labels,
        cellLoc='center',
        loc='center'
    )
    
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    
    # Scale the table to make it bigger (scale factors can be adjusted)
    table.scale(1.5, 1.5)
    
    # Customize the table's appearance: header style and alternating row colors
    for (row, col), cell in table.get_celld().items():
        # Header row formatting
        if row == 0:
            cell.set_text_props(weight='bold', color='white')
            cell.set_facecolor('#40466e')  # dark blue-gray header background
        else:
            # Alternate row colors for readability
            if row % 2 == 0:
                cell.set_facecolor('#f1f1f2')  # light gray
            else:
                cell.set_facecolor('white')
    
    # Save the data if needed
    save_data(summary_df, 'mergersum.dat')
    
    
    
    
    filtered_df = df[df['mergMass1(26)'] > 0.001]
    
    
    
    uniquemergers = sorted(filtered_df['compType'].unique())
    
    
    grouped = filtered_df.groupby(['time[Myr]', 'compType']).size().unstack(fill_value=0)
    
    
    cumulative_counts_df = grouped.cumsum()
    
    
    for comptype in uniquemergers:
        if comptype not in cumulative_counts_df.columns:
            cumulative_counts_df[comptype] = 0
    
    
    cumulative_counts_df = cumulative_counts_df.sort_index()
    x_values = cumulative_counts_df.index
    
 
    for i, comptype in enumerate(cumulative_counts_df.columns):
        y_values = cumulative_counts_df[comptype]
        
    
        if (grouped[comptype] > 0).any():
            last_encounter_time = grouped[comptype][grouped[comptype] > 0].index[-1]
        else:
            last_encounter_time = x_values[0]
        
        
        step = x_values[1] - x_values[0] if len(x_values) > 1 else 1
        
 
        extended_x = np.append(x_values[x_values <= last_encounter_time], last_encounter_time + step)
        extended_y = np.append(y_values[x_values <= last_encounter_time], 0)
        
   
        axs[2, 0].plot(extended_x, extended_y, 
                       label=f'Comp Type {comptype}', 
                       color=color_palette[i % len(color_palette)])
    
 
    axs[2, 0].set_title("Mergers with Different Stellar Types Vs Time")
    axs[2, 0].set_xlabel("Time (Myr)")
    axs[2, 0].set_ylabel("Cumulative Count")
    axs[2, 0].legend(loc='upper left', bbox_to_anchor=(1, 1.5))
    axs[2, 0].grid(True)
    axs[2, 0].set_yscale('log')
    
    plt.tight_layout()
    
    


  
    event_counts_4 = filtered_df3['lineType(1)'].value_counts()
    
 
    palette = sns.color_palette("Set1", n_colors=len(event_counts_4.index))
    
    sns.barplot(
        x=event_counts_4.index,
        y=event_counts_4.values,
        ax=axs[2, 1])
    # sns.barplot(
    #     x=event_counts_4.index,
    #     y=event_counts_4.values,
    #     ax=axs[2, 1],
    #     hue=event_counts_4.index,
    #     palette=palette,
    #     dodge=False,
    #     legend=False
    # )
    
    axs[2, 1].set_title("IMBH-Compact Object Mergers")
    axs[2, 1].set_xlabel("Event Type")
    axs[2, 1].set_ylabel("Frequency")
    axs[2, 1].set_yscale('log')
    
    
    handles = [
        mpatches.Patch(color=palette[i], label=str(count))
        for i, count in enumerate(event_counts_4.values)
    ]
    
    # axs[2, 1].legend(handles=handles, loc='upper left', bbox_to_anchor=(1, 1))

    plt.suptitle("IMBH MERGER PROFILE", fontsize=16)
    
    plt.tight_layout()
    
    plt.savefig('subplots2.png', dpi=300)
    
    

    fig, ax = plt.subplots(figsize=(10, 6))
    
   
    filtered_df2 = filtered_df2.copy()
    
  
    filtered_df2['MassChange'] = filtered_df2['massNew[Msun]'] - filtered_df2['massOld[Msun](10)']
    mass_increase_points2 = filtered_df2[filtered_df2['MassChange'] !=0].copy()
    
  
    mass_increase_points2['stellar_type'] = mass_increase_points2['compType'].map(config.k_values)
  
    filtered_df2['stellar_type'] = filtered_df2['compType'].map(config.k_values)
    
   
    ax.set_facecolor('white')
    ax.grid(True, color='white', linestyle='-', linewidth=0.5)
    
   
    ax.plot(filtered_df2['time[Myr]'], filtered_df2['massNew[Msun]'], color='k')
    

    color_palette = [
        'crimson',       # #DC143C
        'royalblue',     # #4169E1
        'goldenrod',     # #DAA520
        'mediumseagreen',# #3CB371
        'tomato',        # #FF6347
        'slateblue',     # #6A5ACD
        'darkorange',    # #FF8C00
        'mediumorchid',  # #BA55D3
        'lightcoral',    # #F08080
        'seagreen',      # #2E8B57
        'slategray',     # #708090
        'peru',          # #CD853F
        'steelblue',     # #4682B4
        'sandybrown'     # #F4A460
    ]
    
   
    unique_stellar_types = filtered_df2['stellar_type'].unique()
    color_dict = {stellar_type: color_palette[i % len(color_palette)]
                  for i, stellar_type in enumerate(unique_stellar_types)}
 
    sns.scatterplot(data=mass_increase_points2,
                    x='time[Myr]', y='massNew[Msun]',
                    hue='stellar_type',     
                    style='lineType(1)',      
                    palette=color_palette,    
                    s=50, alpha=0.7, zorder=2,
                    legend='full',
                    ax=ax)                 
    
  
    ax.set_title("Mass Evolution in Mergers Over Time with Each Stellar Type")
    ax.set_xlabel("Time (Myr)")
    ax.set_ylabel("Mass (Msun)")
   # ax.set_xlim([0, 2000])
    ax.grid(True,color='k')
    
  
    leg = ax.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=True)
    leg.get_frame().set_facecolor('gray')
    leg.get_frame().set_edgecolor('white')
    plt.setp(leg.get_texts(), color='white')
    
 
    def on_click(event):
        if event.inaxes == ax:
          
            distances = (mass_increase_points2['time[Myr]'] - event.xdata)**2 + \
                        (mass_increase_points2['massNew[Msun]'] - event.ydata)**2
            closest_index = distances.idxmin()
            event_info = mass_increase_points2.loc[closest_index]
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
    plt.tight_layout()
    plt.savefig('mergers.png', dpi=300)
    plt.show()

    
    return
