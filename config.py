import os
import json


def update_json_file(file_path, key, value):
    """
    Updates or adds a key-value pair in a JSON file.

    Parameters:
    file_path (str): The path of the JSON file to be updated.
    key (str): The key to be updated or added.
    value: The value to be set for the key.
    """
    try:
        # Read the existing data
        data = {}
        with open(file_path, 'r') as file:
            data = json.load(file)

        # Update or add the key-value pair
        data[key] = value

        # Write the updated data back to the file
        with open(file_path, 'w') as file:
            json.dump(data, file, indent=4)

        print(f"Updated {file_path} successfully.")

    except FileNotFoundError:
        print(f"The file {file_path} was not found.")
    except json.JSONDecodeError:
        print("Error decoding JSON from the file.")
    except Exception as e:
        print(f"An error occurred: {e}")

def read_json_file(file_path):
    """
    Reads a JSON file and returns its contents.

    Parameters:
    file_path (str): The path of the JSON file to be read.

    Returns:
    dict: The contents of the JSON file.
    """
    try:
        # Open and read the file
        with open(file_path, 'r') as file:
            data = json.load(file)
        
        return data

    except FileNotFoundError:
        print(f"The file {file_path} was not found.")
        return None
    except json.JSONDecodeError:
        print("Error decoding JSON from the file.")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

# Example Usage
# file_path = 'data.json'
# data = read_json_file(file_path)
# if data is not None:
#     print("Data read from the file:", data)

if __name__ == '__main__':
    project_dir = '/'.join(__file__.split('/')[:-1]) + '/'
    #write the global directory to project_directory json variable
    update_json_file(project_dir+"default_files/default_input.json","project_directory",project_dir)
    update_json_file(project_dir+"default_files/default_input.json","data_directory",project_dir[:-1]+"_data/")

    os.system("git submodule init")
    os.system("git submodule update")
    #make project folders
    # if not os.path.exists("jobs/"):
    #     os.mkdir("jobs/")

    