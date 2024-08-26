import pandas as pd
import numpy as np
import config

def get_cumulative_counts(include='both'):
    if include != 'both':
        raise ValueError("This function is designed to work only with 'both' option.")

    data = config.data
    k_values = config.k_values
    cumulative_counts = {stype: 0 for stype in k_values.values()}

    # Handle both single and binary stars
    if 0 in k_values:
        # For stars where ikb == 0 (could be considered single stars in the previous context)
        single_data = data[data['ikb'] == 0]
        
        # Reset index to avoid alignment issues
        single_data = single_data.reset_index(drop=True)

        # Rows where either sm1 or sm2 is zero
        mask = (single_data['sm1'] != 0) | (single_data['sm2'] != 0)
        
        # Set ik1 or ik2 to NaN if corresponding sm1 or sm2 is zero
        single_data.loc[single_data['sm1'] == 0, 'ik1'] = np.nan
        single_data.loc[single_data['sm2'] == 0, 'ik2'] = np.nan
        
        # Apply mask to include only rows where sm1 or sm2 is non-zero
        single_data = single_data[mask].reset_index(drop=True)
        
        # Set sm1 or sm2 to NaN if they are zero
        single_data.loc[single_data['sm1'] == 0, 'sm1'] = np.nan
        single_data.loc[single_data['sm2'] == 0, 'sm2'] = np.nan
        
        single_data['single_mass'] = single_data['sm1'].combine_first(single_data['sm2'])
        single_data = single_data[~single_data['single_mass'].isna()]

        # Update counts for single data
        single_stars_ik1 = single_data[single_data['ik1'].notna()]['ik1'].map(k_values).value_counts()
        single_stars_ik2 = single_data[single_data['ik2'].notna()]['ik2'].map(k_values).value_counts()

        for stype, count in single_stars_ik1.items():
            cumulative_counts[stype] += count
        for stype, count in single_stars_ik2.items():
            cumulative_counts[stype] += count

    # For stars where ikb != 0 (could be considered binary stars in the previous context)
    binary_data = data[data['ikb'] != 0]
    
    # Reset index to avoid alignment issues
    binary_data = binary_data.reset_index(drop=True)

    # Rows where either sm1 or sm2 is zero
    mask = (binary_data['sm1'] != 0) | (binary_data['sm2'] != 0)
    
    # Set ik1 or ik2 to NaN if corresponding sm1 or sm2 is zero
    binary_data.loc[binary_data['sm1'] == 0, 'ik1'] = np.nan
    binary_data.loc[binary_data['sm2'] == 0, 'ik2'] = np.nan
    
    # Apply mask to include only rows where sm1 or sm2 is non-zero
    binary_data = binary_data[mask].reset_index(drop=True)
    
    # Set sm1 or sm2 to NaN if they are zero
    binary_data.loc[binary_data['sm1'] == 0, 'sm1'] = np.nan
    binary_data.loc[binary_data['sm2'] == 0, 'sm2'] = np.nan
    
    binary_data['single_mass'] = binary_data['sm1'].combine_first(binary_data['sm2'])
    binary_data = binary_data[~binary_data['single_mass'].isna()]

    # Update counts for binary data
    binary_stars_ik1 = binary_data[binary_data['ik1'].notna()]['ik1'].map(k_values).value_counts()
    binary_stars_ik2 = binary_data[binary_data['ik2'].notna()]['ik2'].map(k_values).value_counts()

    for stype, count in binary_stars_ik1.items():
        cumulative_counts[stype] += count
    for stype, count in binary_stars_ik2.items():
        cumulative_counts[stype] += count

    # Convert cumulative counts to a DataFrame
    cumulative_df = pd.DataFrame(list(cumulative_counts.items()), columns=['Stellar Type', 'Cumulative Count'])
    
    return cumulative_df
