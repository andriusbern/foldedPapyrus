import os
from urllib import request
import zipfile
import argparse


SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, 'data')
DATASET_URL_LITE = f'https://surfdrive.surf.nl/files/index.php/s/L0BWq8tw3vHSrOC/download?path=%2F'
DATASET_URL_FULL = f'https://surfdrive.surf.nl/files/index.php/s/P8XXnoGonUCMLSr/download?path=%2F'

def download_protein(pid, full_data=False, unzip=True):
    """Downloads the processed protein sequence from online storage"""

    print(f'Downloading protein sequence {pid}...', end='\r')

    url = DATASET_URL_FULL if full_data else DATASET_URL_LITE
    data_url = url + f'{pid}.zip'
    local_file = os.path.join(DATA_DIR, f'{pid}.zip')
    open(local_file, 'a').close()

    try:
        request.urlretrieve(data_url, local_file)
        if unzip:
            unzip_file(local_file)
            
    except:
        print(f'Files for protein {pid} were not found.\n')


def unzip_file(filename, remove_original=True):
    """Unzips a file"""

    with zipfile.ZipFile(filename, 'r') as zip_ref:
        zip_ref.extractall(DATA_DIR)
    if remove_original:
        os.remove(filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pid', type=str, default=None)
    parser.add_argument('-f', '--pid_file', type=str, default=None)
    parser.add_argument('-a', '--all', action='store_true')
    parser.add_argument('-x', '--full', action='store_true')
    args = parser.parse_args()

    print('Downloading foldedPapyrus data...')

    if args.pid:
        download_protein(args.pid, full_data=args.full)
    
    elif args.pid_file:
        with open(args.pid_file, 'r') as f:
            lines = f.readlines()
            print(f'Downloading {len(lines)} proteins listed in {args.pid_file}...\n')
            for line in lines:
                accession = line.strip().split('_')[0]
                download_protein(accession, full_data=args.full)
    
    elif args.all:
        with open('utils/processed_pids.csv', 'r') as f:
            for pid in f.readlines():
                download_protein(pid.strip(), full_data=args.full)
