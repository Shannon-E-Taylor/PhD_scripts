# Measure density, volumes etc

import subprocess


# Takes 10mins to 1hr to run
print('Counting nuclei')
subprocess.run(f'python count_nuclei.py')
# Runs fast (minutes)
print('Measuring volume')
subprocess.run(f'python measure_volumes.py')

print('Extracting cell numbers')
subprocess.run(f'python extract_cell_numbers.py')

# Slow (hours)
# subprocess.run(f'python measure_density_across_psm.py')