import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns

# Ask user for the file path
file_path = input("Please enter the path to the history file: ")

# Read the file and extract lines
with open(file_path, 'r') as file:
    history_content = file.readlines()

# Find the line with column names (line 188 in this case)
column_line = history_content[187].strip()

# Remove the leading `#` symbol and split the rest of the columns
columns = column_line.lstrip('#').split()

# Replace the incorrect column name `#` with `lineType(1)`
if columns[0] == '#':
    columns[0] = 'lineType(1)'

# Read the data into a DataFrame, skipping the metadata lines and using the corrected columns
df = pd.read_csv(file_path, delim_whitespace=True, skiprows=188, comment='#', names=columns)

# Print the head of the DataFrame to check the data
print("DataFrame Head:")
print(df)

print("Column Names:")
print(df.columns.tolist())
# # Ensure 'eOld' and 'eNew' are treated as floats with scientific notation
# df['eOld'] = df['eOld'].astype(float)
# df['eNew'] = df['eNew'].astype(float)



# Identify points where the mass increased
df['MassChange'] = df['massNew[Msun]'] - df['massOld[Msun](10)']
mass_increase_points = df[df['MassChange'] > 0]

# Set up plot style
sns.set(style="whitegrid")

# Create a grid for visualizations (2x2 layout)
fig, axs = plt.subplots(2, 2, figsize=(14, 10))

### Subplot 1: Mass Evolution Over Time with Significant Events

# Plot mass evolution line
axs[0, 0].plot(df['time[Myr]'], df['massNew[Msun]'], label='Mass Evolution', color='b')

# Plot significant mass increase points as solid red points
increase_points = axs[0, 0].scatter(mass_increase_points['time[Myr]'], 
                                    mass_increase_points['massNew[Msun]'], 
                                    color='r', s=1, label='Mass Increase Events', zorder=5)

# Set labels and title
axs[0, 0].set_title("Mass Evolution Over Time with Significant Events")
axs[0, 0].set_xlabel("Time (Myr)")
axs[0, 0].set_ylabel("Mass (Msun)")
axs[0, 0].set_xlim([0, 2000])
# Add legend
axs[0, 0].legend()

def on_click(event):
    if event.inaxes == axs[0, 0]:  # Ensure we're clicking on the correct plot
        # Get the index of the closest point
        closest_index = ((mass_increase_points['time[Myr]'] - event.xdata)**2 + 
                         (mass_increase_points['massNew[Msun]'] - event.ydata)**2).idxmin()

        # Fetch the corresponding row from the DataFrame
        event_info = mass_increase_points.loc[closest_index]

        # Create an elaborated output for this event (displaying all columns)
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
        """
        print(output_text)


# Connect the click event to the function
fig.canvas.mpl_connect('button_press_event', on_click)

### Subplot 2: Event Type Frequency

event_counts = df['lineType(1)'].value_counts()
sns.barplot(x=event_counts.index, y=event_counts.values, ax=axs[0, 1])
axs[0, 1].set_title("Event Type Frequency")
axs[0, 1].set_xlabel("Event Type")
axs[0, 1].set_ylabel("Frequency")
axs[0, 1].tick_params(axis='x', rotation=45)

### Subplot 3: aOld and aNew Over Time

#axs[1, 0].plot(df['time[Myr]'], df['aOld[Rsun]'], label='aOld[Rsun]', color='b')
axs[1, 0].plot(df['time[Myr]'], df['aNew[Rsun]'], label='aNew[Rsun]', color='r')
axs[1, 0].set_title("Semi-major Axis Evolution Over Time")
axs[1, 0].set_xlabel("Time (Myr)")
axs[1, 0].set_ylabel("Semi-major Axis (Rsun)")
axs[1, 0].legend()

### Subplot 4: eOld and eNew Over Time

#axs[1, 1].plot(df['time[Myr]'], df['eOld(16)'], label='eOld', color='b')
axs[1, 1].plot(df['time[Myr]'], df['eNew'], label='eNew', color='r')
axs[1, 1].set_title("Orbital Eccentricity Evolution Over Time")
axs[1, 1].set_xlabel("Time (Myr)")
axs[1, 1].set_ylabel("Orbital Eccentricity")
axs[1, 1].legend()

# Adjust layout and show the plots
plt.tight_layout()
plt.show()
