import os
import pandas as pd
import matplotlib.pyplot as plt

# Dictionary of y-column names with keys starting from 0
y_columns = {
    0: 'Low-mass main-sequence star',
    1: 'Main-sequence star',
    2: 'Hertzsprung-gap star',
    3: 'First giant branch star',
    4: 'Core helium burning star',
    5: 'Early asymptotic giant branch star',
    6: 'Thermally pulsing asymptotic giant branch star',
    7: 'Naked helium star MS',
    8: 'Naked helium star Hertzsprung gap',
    9: 'Naked helium star giant branch',
    10: 'Helium white dwarf',
    11: 'Carbon-oxygen white dwarf',
    12: 'Oxygen-neon white dwarf',
    13: 'Neutron star',
    14: 'Black hole',
    15: 'Massless supernova'
}

def load_data(file_path):
    """Load data from a .dat file with tab delimiter and extract the half-mass radius."""
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        
        # Extract half-mass radius from the 4th line
        half_mass_radius = None
        for line in lines:
            if line.startswith('# Half-Mass Radius:'):
                half_mass_radius = float(line.split(':')[1].split()[0])
                break
        
        if half_mass_radius is None:
            raise ValueError(f"Half-Mass Radius not found in file: {file_path}")

        # Load the data with pandas, specifying tab delimiter and skipping metadata lines
        data = pd.read_csv(file_path, sep='\t', comment='#')

        # Print columns for debugging
        print(f"Columns in {os.path.basename(file_path)}: {data.columns.tolist()}")
        
        return data, half_mass_radius

    except pd.errors.ParserError as e:
        print(f"Error parsing file {file_path}: {e}")
    except Exception as e:
        print(f"An error occurred while loading file {file_path}: {e}")
    return None, None

def plot_cumulative_vs_r(data_files, x_axis_column, selected_y_columns):
    """Plot cumulative counts vs. specified x-axis column for each data file."""
    if x_axis_column not in ['r', 'r/rh']:
        raise ValueError("x_axis_column must be 'r' or 'r/rh'")

    num_files = len(data_files)
    num_cols = 2  # Number of columns in the plot grid
    num_rows = (num_files + num_cols - 1) // num_cols  # Calculate number of rows needed

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(12, 5 * num_rows))
    axes = axes.flatten()  # Flatten the 2D array of axes

    for i, file_path in enumerate(data_files):
        file_name = os.path.basename(file_path)
        data, half_mass_radius = load_data(file_path)

        if data is None:
            print(f"Skipping file {file_name} due to loading issues.")
            continue

        # Check if the x-axis column exists in the data
        if x_axis_column not in data.columns:
            print(f"Column {x_axis_column} not found in file: {file_name}")
            continue

        ax = axes[i]

        # Plot each of the selected y_columns that are present in the data
        for y_key in selected_y_columns:
            y_column = y_columns[y_key]
            if y_column in data.columns:
                if pd.api.types.is_numeric_dtype(data[x_axis_column]) and pd.api.types.is_numeric_dtype(data[y_column]):
                    ax.plot(data[x_axis_column], data[y_column], label=y_column)
                else:
                    print(f"Column {y_column} or {x_axis_column} in file {file_name} is not numeric.")
            else:
                print(f"Column {y_column} not found in file {file_name}.")
        
        ax.set_yscale('log')
        ax.set_xscale('log')
        # Add vertical line at the half-mass radius
        ax.axvline(x=half_mass_radius, color='red', linestyle='--', label=f'Half-Mass Radius')
        
        # Annotate the half-mass radius value on the plot
        ymin, ymax = ax.get_ylim()
        ax.text(half_mass_radius, ymin + (ymax - ymin) * 0.05, f'{half_mass_radius:.2f} pc', color='red', ha='right')

        ax.set_xlabel(x_axis_column)
        ax.set_ylabel('Cumulative Counts')
        ax.set_title(f'Cumulative Counts vs. {x_axis_column} for {file_name}')
        ax.grid(True)

    # Create a single consolidated legend outside the plots
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, title='Stellar Types', loc='upper right', bbox_to_anchor=(1.0, 1.0))

    # Hide any unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.show()

if __name__ == "__main__":
    # Display the y_column options to the user
    print("Available stellar types for plotting:")
    for key, value in y_columns.items():
        print(f"{key}: {value}")

    # User input for selecting which stellar types to plot
    selected_columns_input = input("Enter the numbers corresponding to the stellar types you want to plot, separated by commas (press Enter to select all): ")
    if selected_columns_input.strip() == "":
        # If no input is provided, select all stellar types
        selected_y_columns = list(y_columns.keys())
    else:
        selected_y_columns = [int(x.strip()) for x in selected_columns_input.split(',')]

    # User input for file paths and x-axis column
    file_input = input("Enter the paths to data files, separated by commas: ")
    data_files = [file.strip() for file in file_input.split(',')]
    x_axis_column = input("Enter the column name for the x-axis ('r' or 'r/rh'): ")
    
    if x_axis_column not in ['r', 'r/rh']:
        print("Invalid x-axis column. Please enter 'r' or 'r/rh'.")
    else:
        plot_cumulative_vs_r(data_files, x_axis_column, selected_y_columns)
