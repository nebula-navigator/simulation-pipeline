import os

def validate_file_path(file_path):
    """
    Validate a given file path

    Args:
        file_path (str): The path of the file to validate.    
    Raises:
        ValueError: If the file path is invalid, not accessible, or does not exist.
    """
    if not file_path:
        raise ValueError("Error: The file path cannot be empty.")
    if not isinstance(file_path, str):
        raise ValueError("Error: The file path must be a string.")
    if any(char in file_path for char in ['<', '>', ':', '"', '/', '\\', '|', '?', '*']):
        raise ValueError("Error: The file path contains invalid characters.")
    if not os.path.isfile(file_path):
        raise ValueError(f"Error: The file '{file_path}' does not exist.")
    if not os.access(file_path, os.R_OK):
        raise ValueError(f"Error: The file '{file_path}' is not accessible or is not readable.")
