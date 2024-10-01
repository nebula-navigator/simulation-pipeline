import re
import pandas as pd
import numpy as np

def bhinfrad(bh_mass):
    
    def make_unique(column_names):
        seen = {}
        for i, name in enumerate(column_names):
            if name in seen:
                seen[name] += 1
                column_names[i] = f"{name}_{seen[name]}"
            else:
                seen[name] = 0
        return column_names

    def extract_column_names(header_lines):
        column_names = []
        for line in header_lines:
            match = re.search(r'"name": "([^"]+)"', line)
            if match:
                column_names.append(match.group(1))
        return column_names

    def read_system_dat(file_path):
        
        try:
            with open(file_path, 'r') as file:
                lines = file.readlines()

            header_lines = [line.strip('# ') for line in lines[:2450] if line.startswith('#') and line.strip('# ').strip()]
            column_names = extract_column_names(header_lines)

            column_names = make_unique(column_names)
            data = pd.read_csv(file_path, sep='\s+', skiprows=2451, names=column_names)

            return data, column_names
        except Exception as e:
            print(f"Error in read_system_dat: {e}")
            return pd.DataFrame(), []

    # Constants
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    solar_mass = 1.989e30  # Solar mass in kg
    km_to_m = 1e3  # Conversion from km to m
    parsec_to_m = 3.086e16  # Conversion from parsecs to meters
    
    file_path = input("Enter path to system.dat for calculations : ")
    data, columns = read_system_dat(file_path)

    if data.empty:
        print("No data available to calculate velocity dispersion.")
        return None

    if 'vcn' not in data.columns:
        print("'vcn' column not found in the data.")
        return None
    
    time_period_myr = float(input("Enter the time period in Myr for velocity dispersion calculation: "))
    
    if 'tphys' not in data.columns:
        print("'tphys' (time) column not found in the data.")
        return None

    time_step_index = (np.abs(data['tphys'] - time_period_myr)).idxmin()
    
    lower_bound = max(0, time_step_index - 3)
    upper_bound = min(len(data) - 1, time_step_index + 3)

    time_values = data['tphys'].iloc[lower_bound:upper_bound + 1]
    vcn_values = data['vcn'].iloc[lower_bound:upper_bound + 1]
    print("\nSelected time and vcn values with their indices:")
    for idx, (t_value, vcn_value) in zip(time_values.index, zip(time_values, vcn_values)):
        print(f"Index: {idx}, Time: {t_value:.4f} Myr, vcn: {vcn_value:.4f}")
    print("calculating mean vd...")
    velocity_dispersion = np.mean(vcn_values)  
    print(f"Mean velocity dispersion at {time_period_myr} Myr : ", velocity_dispersion)
    
    bh_mass_kg = bh_mass * solar_mass
    velocity_dispersion_m_s = velocity_dispersion * km_to_m

    
    influence_radius_m = G * bh_mass_kg / (velocity_dispersion_m_s ** 2)

    
    influence_radius = influence_radius_m / parsec_to_m
    
    if influence_radius is not None:
       print(f"The influence radius (velocity dispersion method) for {bh_mass} solar mass black hole is: {influence_radius:.4f} pc")


    return influence_radius

