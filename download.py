import os
from urllib import request
import zipfile
import argparse


SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, 'data')
DATASET_URL = f'https://surfdrive.surf.nl/files/index.php/s/O6uudK28hOfjZxJ/download?path=%2F'


def download_protein(pid, unzip=True):
    """Downloads processed protein sequence from Surfdrive"""

    print(f'Downloading protein sequence {pid}...', end='\r')
    remote_url = os.path.join(DATASET_URL, f'{pid}.zip')
    local_file = os.path.join(DATA_DIR, f'{pid}.zip')
    open(local_file, 'a').close()

    try:
        request.urlretrieve(remote_url, local_file)
        if unzip:
            unzip_file(local_file)
            os.remove(local_file)
    except:
        print(f'Protein sequence {pid} could not be downloaded.')


def unzip_file(filename):
    """Unzips a file"""

    with zipfile.ZipFile(filename, 'r') as zip_ref:
        zip_ref.extractall(DATA_DIR)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pid', type=str, default=None)
    parser.add_argument('-f', '--pid_file', type=str, default=None)
    parser.add_argument('-a', '--all', action='store_true')
    args = parser.parse_args()

    if args.pid:
        download_protein(args.pid)
    
    elif args.pid_file:
        with open(args.pid_file, 'r') as f:
            for line in f:
                download_protein(line.strip())
    
    elif args.all:
        with open('utils/processed_pids.csv', 'r') as f:
            for pid in f.readlines():
                download_protein(pid.strip())
