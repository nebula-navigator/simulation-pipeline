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
        single_data = data[data['ikb'] == 0]
        
        single_data = single_data.reset_index(drop=True)

        mask = (single_data['sm1'] != 0) | (single_data['sm2'] != 0)
    
        single_data.loc[single_data['sm1'] == 0, 'ik1'] = np.nan
        single_data.loc[single_data['sm2'] == 0, 'ik2'] = np.nan
    
        single_data = single_data[mask].reset_index(drop=True)
        

        # Update counts for single data
        single_stars_ik1 = single_data[single_data['ik1'].notna()]['ik1'].map(k_values).value_counts()
        single_stars_ik2 = single_data[single_data['ik2'].notna()]['ik2'].map(k_values).value_counts()

        for stype, count in single_stars_ik1.items():
            cumulative_counts[stype] += count
        #print("low mass main sequence ik1:",cumulative_counts[0])
        for stype, count in single_stars_ik2.items():
            cumulative_counts[stype] += count
        #print("low mass main sequence ik2:",cumulative_counts[0])
    binary_data = data[data['ikb'] != 0]
    
    binary_data = binary_data.reset_index(drop=True)

    mask = (binary_data['sm1'] != 0) | (binary_data['sm2'] != 0)
    
    # Set ik1 or ik2 to NaN if corresponding sm1 or sm2 is zero
    binary_data.loc[binary_data['sm1'] == 0, 'ik1'] = np.nan
    binary_data.loc[binary_data['sm2'] == 0, 'ik2'] = np.nan
    
    

    # Update counts for binary data
    binary_stars_ik1 = binary_data[binary_data['ik1'].notna()]['ik1'].map(k_values).value_counts()
    binary_stars_ik2 = binary_data[binary_data['ik2'].notna()]['ik2'].map(k_values).value_counts()

    for stype, count in binary_stars_ik1.items():
        cumulative_counts[stype] += count
    #print("low mass main sequence ik1 (binary):",cumulative_counts[0])
    for stype, count in binary_stars_ik2.items():
        cumulative_counts[stype] += count
    #print("low mass main sequence ik2 (binary):",cumulative_counts[0])

    # Convert cumulative counts to a DataFrame
    cumulative_df = pd.DataFrame(list(cumulative_counts.items()), columns=['Stellar Type', 'Cumulative Count'])
    
    return cumulative_df
