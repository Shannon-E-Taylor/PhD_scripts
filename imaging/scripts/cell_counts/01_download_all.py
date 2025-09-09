##
# run this code to run all downloads and analyses
import getpass
OMEROPASS = getpass.getpass()

import subprocess

# Download data from OMERO
# requires access to OMERO etc
subprocess.run(f'python download_data.py {OMEROPASS}')
# subprocess.run(f'python download_somite_annotations.py {OMEROPASS}')

