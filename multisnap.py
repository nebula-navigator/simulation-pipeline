import os
import pandas as pd
import matplotlib.pyplot as plt

def load_data(file_path):
    """Load data from a file, skipping metadata lines and using file-provided column names."""
    try:
        # First, find where the metadata ends and the data begins
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Find the index of the first non-comment line (the column headers)
        for idx, line in enumerate(lines):
            if not line.startswith('#'):
                header_idx = idx
                break

        # Load the data with pandas, skipping the metadata lines
        data = pd.read_csv(file_path, delim_whitespace=True, skiprows=header_idx, header=0)
        return data

    except pd.errors.ParserError as e:
        print(f"Error parsing file {file_path}: {e}")
    except Exception as e:
        print(f"An error occurred while loading file {file_path}: {e}")
    return None

def drop_constant_columns(data):
    """Drop Y columns that become constant (not changing anymore)."""
    for col in data.columns:
        if data[col].nunique() == 1:
            data[col] = 0
    return data

def plot_cumulative_vs_r(data_files, x_axis_column):
    """Plot cumulative counts vs. specified x-axis column for each data file in a grid layout."""
    if x_axis_column not in ['r', 'r/rh']:
        raise ValueError("x_axis_column must be 'r' or 'r/rh'")

    # Calculate grid size
    num_files = len(data_files)
    num_cols = 2
    num_rows = (num_files + 1) // num_cols

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, 5 * num_rows))
    axes = axes.flatten()  # Flatten the 2D array of axes

    for i, file_path in enumerate(data_files):
        file_name = os.path.basename(file_path)
        data = load_data(file_path)

        if data is None:
            print(f"Skipping file {file_name} due to loading issues.")
            continue

        # Check if the x-axis column exists in the data
        if x_axis_column not in data.columns:
            print(f"Column {x_axis_column} not found in file: {file_name}")
            continue

        # Drop constant columns
        data = drop_constant_columns(data)

        # Filter out x-axis and non-numeric columns for y-axis
        y_columns = [col for col in data.columns if col != x_axis_column]
        if not y_columns:
            print(f"No y columns found in file: {file_name}")
            continue

        ax = axes[i]
        for y_column in y_columns:
            ax.plot(data[x_axis_column], data[y_column], label=y_column)

        ax.set_xlabel(x_axis_column)
        ax.set_ylabel('Cumulative Counts')
        ax.set_title(f'Cumulative Counts vs. {x_axis_column} for {file_name}')

        # Create a legend for the y columns
        ax.legend(title='Stellar Types', loc='upper right', bbox_to_anchor=(1.0, 1.0))

        ax.grid(True)

    # Hide any unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    file_input = input("Enter the paths to data files, separated by commas: ")
    data_files = [file.strip() for file in file_input.split(',')]
    x_axis_column = input("Enter the column name for the x-axis ('r' or 'r/rh'): ")
    if x_axis_column not in ['r', 'r/rh']:
        print("Invalid x-axis column. Please enter 'r' or 'r/rh'.")
    else:
        plot_cumulative_vs_r(data_files, x_axis_column)
