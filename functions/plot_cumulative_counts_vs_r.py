import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import config
from .calculate_core_radius import calculate_core_radius
from .calculate_half_mass_radius import calculate_half_mass_radius
from .get_luminosity import get_luminosity
import os
import datetime
import numpy as np

def save_data(data, stage_name, metadata=""):
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    file_name = f"{stage_name}_{timestamp}.dat"
    
    
    if not os.path.exists('data_saves'):
        os.makedirs('data_saves')
    
    file_path = os.path.join('data_saves', file_name)
    
    # Preparing the data to be saved with aligned columns
    column_widths = [max(len(str(value)) for value in [col] + data[col].astype(str).tolist()) + 2 for col in data.columns]
    
    with open(file_path, 'w') as f:
        # Writing metadata if provided
        if metadata:
            f.write("# Metadata\n")
            f.write(metadata)
            f.write("\n")
        
        f.write("# Data\n")
        
        # Writing column headers with proper alignment
        for col, width in zip(data.columns, column_widths):
            f.write(f"{col.ljust(width)}")
        f.write('\n')
        
        # Writing data rows with proper alignment
        for _, row in data.iterrows():
            for col, width in zip(data.columns, column_widths):
                f.write(f"{str(row[col]).ljust(width)}")
            f.write('\n')
    
    print(f"Data saved to {file_path}")
    
def plot_cumulative_counts_vs_r(current_data):
    data = current_data
    k_values = config.k_values
    include_options = ['single', 'binary', 'both']
    
    print("Select the type of stars to include:")
    print("1. Single Stars")
    print("2. Binary Stars")
    print("3. Both")
    choice = input("Enter your choice (1-3): ")

    if choice == '1':
        include = 'single'
    elif choice == '2':
        include = 'binary'
    elif choice == '3':
        include = 'both'
    else:
        print("Invalid choice. Defaulting to 'both'.")
        include = 'both'

    print("Available stellar types:")
    for k, stype in k_values.items():
        print(f"{k}: {stype}")

    valid_types_input = input("Enter stellar type numbers to include, separated by commas (e.g., '1,2,3'): ")

    if valid_types_input:
        valid_types = [int(t.strip()) for t in valid_types_input.split(',')]
    else:
        valid_types = list(k_values.keys())
    

    if 'ikb' not in data.columns:
        print("Error: 'ikb' column not found in data, which is required for handling k=0.")
        return

    # Filtering based on the selected include type early on
    if include == 'single':
        filtered_data = data[data['ikb'] == 0].copy()
    elif include == 'binary':
        filtered_data = data[data['ikb'] != 0].copy()
    else:  # both
        filtered_data = data.copy()
    
    print("Data before filtering for k values:")
    print(filtered_data.head()) 
    # metadata = f"# Data before filtering for k values\n# Timestamp: {datetime.datetime.now().isoformat()}\n"
    # save_data(filtered_data, "data_before_filtering", metadata)
    

    # for k=0 (single stars)
    if 0 in valid_types or include == 'single' :
        # rows where either sm1 or sm2 is zero 
        mask = (filtered_data['sm1'] != 0) | (filtered_data['sm2'] != 0)

        # If sm1 is zero, set ik1 to NaN; if sm2 is zero, set ik2 to NaN
        filtered_data.loc[filtered_data['sm1'] == 0, 'ik1'] = np.nan
        filtered_data.loc[filtered_data['sm2'] == 0, 'ik2'] = np.nan

        # Combining columns into a new column 'single_mass' by summing the values of 'sm1' and 'sm2' per row
        filtered_data.loc[:, 'single_mass'] = filtered_data[['sm1', 'sm2']].sum(axis=1)


        # Dropping rows where 'single_mass' is NaN
        filtered_data = filtered_data.dropna(subset=['single_mass'])
       
    
       # mask to only include rows where 'single_mass' is non-NaN
        filtered_data = filtered_data[mask]
    
    if 0 in valid_types or include == 'binary' :
        # rows where either sm1 or sm2 is zero 
        mask = (filtered_data['sm1'] != 0) | (filtered_data['sm2'] != 0)

        # If sm1 is zero, set ik1 to NaN; if sm2 is zero, set ik2 to NaN
        filtered_data.loc[filtered_data['sm1'] == 0, 'ik1'] = np.nan
        filtered_data.loc[filtered_data['sm2'] == 0, 'ik2'] = np.nan

        # Combining columns into a new column 'single_mass' by summing the values of 'sm1' and 'sm2' per row
        filtered_data.loc[:, 'single_mass'] = filtered_data[['sm1', 'sm2']].sum(axis=1)


       # Dropping rows where 'single_mass' is NaN
        filtered_data = filtered_data.dropna(subset=['single_mass'])
      
    
       # masking to only include rows where 'single_mass' is non-NaN
        filtered_data = filtered_data[mask]
        
    if 0 in valid_types or include == 'both' :
        # rows where either sm1 or sm2 is zero 
        mask = (filtered_data['sm1'] != 0) | (filtered_data['sm2'] != 0)

        # If sm1 is zero, set ik1 to NaN; if sm2 is zero, set ik2 to NaN
        filtered_data.loc[filtered_data['sm1'] == 0, 'ik1'] = np.nan
        filtered_data.loc[filtered_data['sm2'] == 0, 'ik2'] = np.nan
        
    
        
        # Combine columns into a new column 'single_mass' by summing the values of 'sm1' and 'sm2' per row
        filtered_data.loc[:, 'single_mass'] = filtered_data[['sm1', 'sm2']].sum(axis=1)


         # Drop rows where 'single_mass' is NaN
        filtered_data = filtered_data.dropna(subset=['single_mass'])
       
    
        
    print("Data after handling special case for k=0 (if applicable):")
    print(filtered_data.head()) 
    # Save data after handling special case for k=0
    #metadata = f"# Data after handling special case for k=0\n# Timestamp: {datetime.datetime.now().isoformat()}\n"
    #save_data(filtered_data, "data_after_handling_k0", metadata)
    
    filtered_data = filtered_data[
        (filtered_data['ik1'].isin(valid_types)) | (filtered_data['ik2'].isin(valid_types))
    ]

    print("Data after filtering for selected k values:")
    print(filtered_data.head()) 
    #Save data after filtering for k values
    # metadata = f"# Data after filtering for selected k values\n# Timestamp: {datetime.datetime.now().isoformat()}\n"
    # save_data(filtered_data, "data_after_filtering", metadata)
    
    #half mass radius
    if include == 'single':
        filtered_data['combined_mass'] = filtered_data['sm1'] + filtered_data['sm2']
        #half_mass_radius = calculate_half_mass_radius(filtered_data, 'combined_mass')
        half_mass_radius = calculate_half_mass_radius(filtered_data, 'single_mass')
    elif include == 'binary':
        filtered_data['combined_mass'] = filtered_data['sm1'] + filtered_data['sm2']
        #half_mass_radius = calculate_half_mass_radius(filtered_data, 'combined_mass')
        half_mass_radius = calculate_half_mass_radius(filtered_data, 'single_mass')
    else:  # both
       #filtered_data['combined_mass'] = filtered_data['sm1'] + filtered_data['sm2']
       #half_mass_radius = calculate_half_mass_radius(filtered_data, 'combined_mass')
        half_mass_radius = calculate_half_mass_radius(filtered_data, 'single_mass')

    filtered_data['effective_lum'] = filtered_data.apply(
        lambda row: get_luminosity(row['sm1'], row['sm2'], row['lum1'], row['lum2']),
        axis=1
    )
    filtered_data['core_radius'] = calculate_core_radius(filtered_data['effective_lum'])
    core_radius_mean = filtered_data['core_radius'].mean()

    print(f"Half-Mass Radius: {half_mass_radius:.2f} pc")
    print(f"Mean Core Radius: {core_radius_mean:.2f} pc")

    # additional filtering based on radial position
    min_radial = input("Enter minimum radial position (pc) or press Enter to skip: ")
    max_radial = input("Enter maximum radial position (pc) or press Enter to skip: ")

    if min_radial:
        min_radial = float(min_radial)
        filtered_data = filtered_data[filtered_data['r'] >= min_radial]

    if max_radial:
        max_radial = float(max_radial)
        filtered_data = filtered_data[filtered_data['r'] <= max_radial]

    if filtered_data.empty:
        print("No data available for the selected filters.")
        return

    x_axis_options = ['r'] + [col for col in filtered_data.columns if col.startswith('r/rh')]
    print("Select the x-axis column:")
    for i, option in enumerate(x_axis_options, 1):
        print(f"{i}. {option}")

    x_column_choice = int(input(f"Enter your choice (1-{len(x_axis_options)}): "))
    x_column = x_axis_options[x_column_choice - 1]
    
    # Combine ik1 and ik2 into a single column 'combined_ik'
    filtered_data['combined_ik'] = filtered_data[['ik1', 'ik2']].apply(lambda x: [val for val in x if pd.notna(val)], axis=1)

    # Explode the list column 'combined_ik' to have one row per stellar type
    filtered_data = filtered_data.explode('combined_ik')

    # Convert the combined_ik column to an integer for easier filtering later
    filtered_data['combined_ik'] = filtered_data['combined_ik'].astype(int)

    cumulative_counts_df = pd.DataFrame()

    # Group by radial position and combined stellar types, and count occurrences
    grouped = filtered_data.groupby([x_column, 'combined_ik']).size().unstack(fill_value=0)

# Cumulative sum for each stellar type across all radial positions
    cumulative_counts_df = grouped.cumsum()

# Filter to only include valid types
    cumulative_counts_df = cumulative_counts_df.loc[:, cumulative_counts_df.columns.isin(valid_types)]

# Fill missing stellar types with 0 cumulative count if they don't exist in the data
    for stype in valid_types:
       if stype not in cumulative_counts_df.columns:
        cumulative_counts_df[stype] = 0

# Rename columns to stellar type names
    cumulative_counts_df.columns = [k_values[stype] for stype in cumulative_counts_df.columns]

# Sort by radial position for proper cumulative plotting
    cumulative_counts_df = cumulative_counts_df.sort_index()

    x_values = cumulative_counts_df.index
    
    plt.figure(figsize=(10, 6))

# Plot each stellar type's cumulative count as a separate line
    for column in cumulative_counts_df.columns:
       plt.plot(x_values, cumulative_counts_df[column], label=f'{column}')

    plt.axvline(x=half_mass_radius, color='r', linestyle='--', label=f'Half-Mass Radius: {half_mass_radius:.2f} pc')
    plt.axvline(x=core_radius_mean, color='b', linestyle='--', label=f'Mean Core Radius: {core_radius_mean:.2f} pc')

    plt.xlabel(f'{x_column} (pc)')
    plt.ylabel('Cumulative Count')
    plt.title('Cumulative Counts vs. Radial Position')
    plt.grid(True)
    plt.legend(title='Stellar Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()  

    plt.show()

    save_data_prompt = input("Would you like to save the filtered data? (y/n): ").strip().lower()

    if save_data_prompt == 'y':
        file_name = input("Enter the file name to save the data (without extension): ").strip() + '.dat'
        metadata = {
            'Selected star inclusion': include,
            'Selected stellar types': valid_types,
            'Half-Mass Radius': f"{half_mass_radius:.2f} pc",
            'Mean Core Radius': f"{core_radius_mean:.2f} pc",
            'Minimum radial position': min_radial,
            'Maximum radial position': max_radial,
            'x_column': x_column,
            'Timestamp': datetime.datetime.now().isoformat()
        }

        if filtered_data.index.duplicated().any():
            print("Warning: Duplicate index values found. Resetting index.")
            filtered_data = filtered_data.reset_index(drop=True)

        print("Calculating r/rh...")
        filtered_data['r/rh'] = filtered_data['r'] / half_mass_radius
        if include == 'single':
            filtered_data['r/rh(single)'] = filtered_data['r/rh']
        elif include == 'binary':
            filtered_data['r/rh(binary)'] = filtered_data['r/rh']
        else:
            filtered_data['r/rh(both)'] = filtered_data['r/rh']

        filtered_data_to_save = filtered_data[[x_column, 'r/rh']].merge(cumulative_counts_df, left_on=x_column, right_index=True)

       
        with open(file_name, 'w') as f:
            # Write metadata
            f.write("# Metadata\n")
            for key, value in metadata.items():
                f.write(f"# {key}: {value}\n")

            # Write data
            filtered_data_to_save.to_csv(f, index=False, sep='\t') 

        print(f"Data saved to {file_name}")
    else:
        print("Data not saved.")

    del filtered_data