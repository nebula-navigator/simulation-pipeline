import os 
def extract_file_name(file_path):
    """
    Extract the file name from the given file path.
    
    Args:
        file_path (str): The path of the file.
        
    Returns:
        str: The name of the file, or an error message if the input is invalid.
    """
    if not file_path:
        raise ValueError("Error: The file path cannot be empty.")
    
    file_name = os.path.basename(file_path)
    
    return file_name

def save_data(df, filename, history_basename):
    if df.empty:
        print("df is empty, skipping save")
        return

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