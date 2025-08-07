import pandas as pd
import numpy as np
import os
import ase
from ase.io import read, write
import argparse

def process_xdatcar(xdatcar_file, time_points, time_step, output_dir):
    """
    Reads a XDATCAR file, and for each time in time_points,
    it finds the corresponding structure and writes it to a CONTCAR file.
    Args:
        xdatcar_file (str): Path to the XDATCAR file.
        time_points (list): List of time steps to extract.
        time_step (int): The time step interval used in the XDATCAR file.
        output_dir (str): The directory to save the output CONTCAR files.
    """
    files_written = 0
    try:
        # Read all structures from the XDATCAR file
        all_structures = read(xdatcar_file, index=':')
        
        # The XDATCAR file contains one structure per time step.
        # We assume a 1-to-1 mapping between index and time.
        for t in time_points:
            structure_index = int(t) # time_steps are 1-based from VASP
            if 0 < structure_index <= len(all_structures):
                structure = all_structures[structure_index - 1]
                output_filename = os.path.join(output_dir, f'CONTCAR_{int(t*time_step)}')
                print(f"Writing structure for time step {int(t*time_step)} to {output_filename}")
                write(output_filename, structure, format='vasp')
                files_written += 1
            else:
                print(f"Error: Time step {int(t*time_step)} is out of bounds. Max steps: {len(all_structures)*time_step}")
    except FileNotFoundError:
        print(f"Error: XDATCAR file not found at '{xdatcar_file}'")
    except Exception as e:
        print(f"An error occurred: {e}")
    print(f"Processing completed. Total files written: {files_written}")

def main():
    parser = argparse.ArgumentParser(description='Process AIMD time-frames based on specified variable in XDATCAR, outputs CONTCAR files for selected time points.')
    parser.add_argument('--sort_variable', type=str, default='E0', help='Variable to sort frames by, default is E0.')
    parser.add_argument('--excel_file', required=True, help='Path to the Excel file with energy data.')
    parser.add_argument('--xdatcar_file', required=True, help='Path to the XDATCAR file.')
    parser.add_argument('--time_step', type=int, default=50, help='Time step used in XDATCAR file, default is 50ps.')
    parser.add_argument('--num_time_points', type=int, default=10, help='Number of time points to select, default is 10 points with lowest energy.')
    parser.add_argument('--min_time', type=int, default=1000, help='Minimum time to start sorting frames, default is 1000ps.')
    parser.add_argument('--output_dir', default='.', help='Directory to save CONTCAR files.')

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Read and process data
    df = pd.read_excel(args.excel_file)
    
    if args.min_time < 0 or args.time_step <= 0 or args.num_time_points <= 0:
        print("Error: All parameters must be positive integers.")
        return

    if args.min_time > max(df['time']):
        print(f"Error: min_time ({args.min_time}) is greater than the maximum time in the data ({max(df['time'])}).")
        return
    
    if args.sort_variable not in df.columns:
        print(f"Error: sort_variable ({args.sort_variable}) is not a valid column in the data. Please check the Excel file.")
        return
    
    if 'time' not in df.columns:
        print(f"Error: 'time' column is not present in the data. Please check the format of Excel file.")
        return

    # Filter dataframe
    df_filtered = df[(df['time'] % args.time_step == 0) & (df['time'] > args.min_time)].sort_values(by=args.sort_variable)

    # Check if num_time_points is valid for the filtered data
    if args.num_time_points > len(df_filtered):
        print(f"Error: num_time_points ({args.num_time_points}) is larger than the number of available frames after filtering ({len(df_filtered)}).")
        return


    # Select time points
    time_selected = df_filtered['time'].iloc[:args.num_time_points].tolist()
    
    # Divide the selected time by time step
    time_steps_for_extraction = [t / args.time_step for t in time_selected]

    # Process XDATCAR
    process_xdatcar(args.xdatcar_file, time_steps_for_extraction, args.time_step, args.output_dir)

if __name__ == '__main__':
    main()
