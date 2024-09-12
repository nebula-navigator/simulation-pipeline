import matplotlib.pyplot as plt
import pandas as pd
import config

def blue_straggler(current_data):
    
    data = current_data
    k_values = config.k_values

    try:
        turn_off_mass = float(input("Enter the main sequence turn-off mass for Blue Stragglers (in solar masses): "))
    except ValueError:
        print("Invalid input. Please enter a numerical value.")
        return


    blue_stragglers = data[(data['ik1'] == 1) | (data['ik2'] == 1)].copy()

    
    blue_stragglers['is_blue_straggler'] = False
    blue_stragglers.loc[blue_stragglers['ik1'] == 1, 'is_blue_straggler'] = blue_stragglers['sm1'] >= turn_off_mass
    blue_stragglers.loc[blue_stragglers['ik2'] == 1, 'is_blue_straggler'] = blue_stragglers['sm2'] >= turn_off_mass

    
    blue_stragglers = blue_stragglers[blue_stragglers['is_blue_straggler']]

   
    if 'r' not in blue_stragglers.columns:
        print("Error: 'r' column not found in the data.")
        return

    
    blue_stragglers_sorted = blue_stragglers.sort_values(by='r')

    
    blue_stragglers_sorted['single_mass'] = blue_stragglers_sorted[['sm1', 'sm2']].sum(axis=1)

    
    blue_stragglers_sorted = blue_stragglers_sorted.dropna(subset=['single_mass'])

  
    blue_straggler_counts = blue_stragglers_sorted.groupby('r').size().cumsum()


    blue_straggler_counts_df = blue_straggler_counts.to_frame(name='Blue Stragglers').reset_index()


    total_blue_stragglers = blue_straggler_counts_df['Blue Stragglers'].iloc[-1]

 
    plt.figure(figsize=(12, 8))
    plt.plot(blue_straggler_counts_df['r'], blue_straggler_counts_df['Blue Stragglers'], label='Blue Stragglers', color='blue')

   
    total_counts_text = f"Total Counts: {total_blue_stragglers}"
    plt.gca().text(0.05, 0.95, total_counts_text, transform=plt.gca().transAxes, 
                   fontsize=12, verticalalignment='top', horizontalalignment='left',
                   bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))

    plt.xlabel('Radial Position (pc)')
    plt.ylabel('Cumulative Count')
    plt.title('Cumulative Counts of Blue Stragglers vs. Radial Position')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

