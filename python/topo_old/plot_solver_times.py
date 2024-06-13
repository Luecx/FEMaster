import os
import re
import matplotlib.pyplot as plt

def extract_values_from_file(filename):
    with open(filename, 'r') as file:
        content = file.read()

        # Use regex to extract the required values
        N = re.search(r'N=(\d+)', content)
        nnz = re.search(r'nnz=(\d+)', content)
        elapsed_time = re.search(r'Elapsed time: (\d+)', content)

        if N and nnz and elapsed_time:
            return {
                'N': int(N.group(1)),
                'nnz': int(nnz.group(1)),
                'Elapsed time': int(elapsed_time.group(1))
            }
        else:
            return None

def extract_and_plot(directory):
    records = []

    # Loop through all files in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".log"):
            path = os.path.join(directory, filename)
            record = extract_values_from_file(path)
            if record:
                records.append(record)

    # Sort by 'nnz'
    records = sorted(records, key=lambda x: x['nnz'])

    # Group by 'nnz' and plot 'Elapsed time'
    grouped = {}
    for record in records:
        if record['nnz'] in grouped:
            grouped[record['nnz']].append(record['Elapsed time'])
        else:
            grouped[record['nnz']] = [record['Elapsed time']]

    for nnz, elapsed_times in grouped.items():
        plt.plot(elapsed_times, label=f'nnz={nnz}')

    plt.xlabel('Files')
    plt.ylabel('Elapsed time')
    plt.legend()
    plt.show()

# Use the function on your directory
extract_and_plot('log/')
