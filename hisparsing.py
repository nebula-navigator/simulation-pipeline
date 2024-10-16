import pandas as pd

# File path
file_path = 'history-1893666.dat'  # Replace with your actual file path

# Read the file and extract lines
with open(file_path, 'r') as file:
    history_content = file.readlines()

# Find the line with column names (line 188 in this case)
column_line = history_content[187].strip()

# Remove the leading `#` symbol and split the rest of the columns
columns = column_line.lstrip('#').split()

# Read the data into a DataFrame, skipping the metadata lines and using the corrected columns
df = pd.read_csv(file_path, delim_whitespace=True, skiprows=188, comment='#', names=columns)

# Display the first few rows to verify
print(df.head())
