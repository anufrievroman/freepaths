"""Example script to run several simulations in which temperature is changed"""

import subprocess

FILENAME = "point_line.py"


def update_config_file(config_file, temperature):
    with open(config_file, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.startswith("OUTPUT_FOLDER_NAME"):
            lines[i] = f"OUTPUT_FOLDER_NAME = {temperature}\n"
        if line.startswith("T "):
            lines[i] = f"T = {temperature}\n"

    with open(config_file, 'w') as f:
        f.writelines(lines)


def main():

    # Define your list of temperatures here:
    temperatures = [100, 200, 300]

    for temp in temperatures:
        update_config_file(FILENAME, temp)
        subprocess.run(['python', '-m', 'freepaths', FILENAME])


if __name__ == "__main__":
    main()
