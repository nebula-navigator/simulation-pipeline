import os
import pandas as pd
import matplotlib.pyplot as plt

# List of y-column names provided by you
y_columns = [
    'Low-mass main-sequence star',
    'Main-sequence star',
    'Hertzsprung-gap star',
    'First giant branch star',
    'Core helium burning star',
    'Early asymptotic giant branch star',
    'Thermally pulsing asymptotic giant branch star',
    'Naked helium star MS',
    'Naked helium star Hertzsprung gap',
    'Naked helium star giant branch',
    'Helium white dwarf',
    'Carbon-oxygen white dwarf',
    'Oxygen-neon white dwarf',
    'Neutron star',
    'Black hole',
    'Massless supernova'
]

def load_data(file_path):
    """Load data from a .dat file with tab delimiter and extract the half-mass radius."""
    try:
        # Read the first few lines to get metadata including the half-mass radius
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

def plot_cumulative_vs_r(data_files, x_axis_column):
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

        # Determine which y_columns are present in the data
        present_y_columns = [col for col in y_columns if col in data.columns]

        if not present_y_columns:
            print(f"No y_columns found in file: {file_name}")
            continue

        # Plot each of the present y_columns
        for y_column in present_y_columns:
            if pd.api.types.is_numeric_dtype(data[x_axis_column]) and pd.api.types.is_numeric_dtype(data[y_column]):
                ax.plot(data[x_axis_column], data[y_column], label=y_column)
            else:
                print(f"Column {y_column} or {x_axis_column} in file {file_name} is not numeric.")

        # Add vertical line at the half-mass radius
        ax.axvline(x=half_mass_radius, color='red', linestyle='--', label=f'Half-Mass Radius: {half_mass_radius:.2f} pc')
        
        # Annotate the half-mass radius value on the plot
        ymin, ymax = ax.get_ylim()
        ax.text(half_mass_radius, ymin + (ymax - ymin) * 0.05, f'{half_mass_radius:.2f} pc', color='red', ha='right')

        ax.set_xlabel(x_axis_column)
        ax.set_ylabel('Cumulative Counts')
        ax.set_title(f'Cumulative Counts vs. {x_axis_column} for {file_name}')
        ax.grid(True)

    # Create a single consolidated legend outside the plots
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, title='Stellar Types', loc='upper right', bbox_to_anchor=(1.1, 1.0))

    # Hide any unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    # Adjust layout to make room for the legend
    plt.tight_layout(rect=[0, 0, 0, 1])  # Adjust the right side to make room for the legend
    plt.show()

if __name__ == "__main__":
    file_input = input("Enter the paths to data files, separated by commas: ")
    data_files = [file.strip() for file in file_input.split(',')]
    x_axis_column = input("Enter the column name for the x-axis ('r' or 'r/rh'): ")
    if x_axis_column not in ['r', 'r/rh']:
        print("Invalid x-axis column. Please enter 'r' or 'r/rh'.")
    else:
        plot_cumulative_vs_r(data_files, x_axis_column)
