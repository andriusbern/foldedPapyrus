import os
import shutil
import torch
import pickle
import numpy as np
from tqdm import tqdm
import pandas as pd
from utils import VocSmiles

LENGTH_PAPYRUS = 61085165


class ProteinDataset:
    """Class for loading and storing AlphaFold embeddings in torch.Tensor format"""

    def __init__(self, data_dir, voc) -> None:
        self.data_dir = data_dir
        self.protein_dir = os.path.join(data_dir, 'proteins')
        self.max_len = 450
        
        # Embedding dict accessible via protein_id
        self.protein_embeddings = self.load_protein_embeddings()
        self.voc = voc
        
    def get_protein_embedding(self, protein_id, embedding_type='single'):
        """
        Returns a torch tensor of the target protein 
        """
        protein_dir = os.path.join(self.data_dir, 'proteins', protein_id)
        embedding_file = os.path.join(protein_dir, f'{embedding_type}.npy')

        np_file = np.load(embedding_file)
        embedding = torch.from_numpy(np_file)
        n, f = embedding.shape
        output = torch.zeros(self.max_len, f) 
        output[:n, :f] = embedding
        return output

    def load_protein_embeddings(self, embedding_type='single'):

        protein_dict = {}
        protein_ids =  [folder for folder in os.listdir(self.protein_dir) 
                        if os.path.isdir(os.path.join(self.protein_dir, folder))]

        for pid in protein_ids:
            protein_dict[pid] = self.get_protein_embedding(pid)
        return protein_dict

        
class ProteinSmilesDataset(torch.utils.data.Dataset, ProteinDataset):
    """Class for storing both AF and papyrus data"""

    def __init__(self, data_dir, dataset_prefix, batch_size=32, **kwargs):
        
        super(ProteinSmilesDataset, self).__init__(data_dir=data_dir, **kwargs)
        self.data_dir = data_dir
        self.dataset_path = os.path.join(data_dir, f'{dataset_prefix}.tsv')
        self.tsv_dataset = self.read_dataset()
        print(f'Dataset: len: {len(self)}')

    def __len__(self):
        return len(self.tsv_dataset)

    def __getitem__(self, idx):
        """
        Returns torch tensors for:
            1) Protein embeddings
            2) Standardized SMILES
            3) Pchembl value scalar
        """

        data = self.tsv_dataset[idx]
        pid, pchembl, tokens = data.split('\t')

        pchembl = torch.tensor(float(pchembl))
        encoded_smiles = self.voc.encode([tokens.strip('\n').split(' ')])
        protein_embeddings = self.protein_embeddings[pid]
        return protein_embeddings, encoded_smiles, pchembl

    def read_dataset(self):
        """
        Reads the dataset, makes sure that embeddings for each protein are available
        """
        
        filtered = []
        with open(self.dataset_path, 'r') as data:
            for line in data.readlines():
                
                pid, pchembl, tokens = line.split('\t')
                if self.protein_embeddings.get(pid) is None:
                    continue
                filtered.append(line)
        
        return filtered 


def clean_dataset(af_output_dir, output_dir, protein_ids=None, zip=False):
    """
    Cleaning
    """

    if protein_ids is None:

        protein_ids = [folder for folder in os.listdir(af_output_dir) \
                       if os.path.isdir(os.path.join(af_output_dir, folder))]

    print(protein_ids)
    print(output_dir)

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'proteins'), exist_ok=True)

    # Process AlphaFold's output
    for i, protein in enumerate(protein_ids):
        try:
            print(f'Processing {protein} | {i}/{len(protein_ids)}', end='\r')
            
            src_dir = os.path.join(af_output_dir, protein)
            protein_dir = os.path.join(output_dir, 'proteins', protein)

            path = os.path.join(src_dir, 'result_model_1_pred_0.pkl')
            pkl = open(path, 'rb')
            result = pickle.load(pkl)
            reps = result['representations']
            single, structure = reps['single'], reps['structure_module']
            pair, msa = reps['pair'], reps['msa']
            plddt = result['plddt']                          # [L]
            distogram_logits = result['distogram']['logits'] # [L x L x 64]
            confidence = result['ranking_confidence']        # [1]
            
            os.makedirs(protein_dir, exist_ok=True)
            np.save(os.path.join(protein_dir, 'single.npy'), single)
            np.save(os.path.join(protein_dir, 'structure.npy'), structure)
            np.save(os.path.join(protein_dir, 'pair.npy'), pair)
            np.save(os.path.join(protein_dir, 'msa.npy'), msa)
            np.save(os.path.join(protein_dir, 'plddt.npy'), plddt)
            np.save(os.path.join(protein_dir, 'distogram_logits.npy'), distogram_logits)
            np.save(os.path.join(protein_dir, 'confidence.npy'), confidence)

            files = ['ranked_0.pdb',
                     'unrelaxed_model_1_pred_0.pdb', 
                     f'{protein}.fasta',
                     'timings.json']

            names = [f'{protein}.pdb',
                     f'{protein}_unrelaxed.pdb',
                     f'{protein}.fasta',
                     'timings.json']
    
            for f, n in zip(files, names):
                shutil.copyfile(os.path.join(src_dir, f), os.path.join(protein_dir, n))

            # If needed zips the output folder and deletes the original    
            if zip:
                shutil.make_archive(output_dir, 'zip', protein_dir)
                shutil.rmtree(protein_dir)

        except:
            pass

def build_dataset(af_output_dir, output_dir, papyrus_file, output_dataset_prefix, 
                  protein_targets_file=None, save_voc=True, process=False):
    """
    For building a compact dataset with only really required data
    1. List all processed proteins
    2. Process the output if needed, move to target dir
    3. Read papyrus line by line

        - if protein in entry is processed, store this as a line
        - at the end each protein has a list of smiles
        - they can be processed via the usual pipeline
    4. Conjoin all the data as
        pid  smiles  pX

    """
    processed_proteins = [folder for folder in os.listdir(af_output_dir) 
                   if os.path.isdir(os.path.join(af_output_dir, folder))]

    if protein_targets_file:
        with open(protein_targets_file, 'r') as prot_list:
            target_proteins = prot_list.readlines()
            target_proteins = [t.strip('\n') for t in target_proteins]
    else:
        target_proteins = processed_proteins
    
    if process:
        clean_dataset(af_output_dir, output_dir, protein_ids=target_proteins)

    # Rescan
    t_dir = os.path.join(output_dir, 'proteins')
    print(t_dir)
    processed_proteins = [folder for folder in os.listdir(t_dir) 
                          if os.path.isdir(os.path.join(t_dir, folder))]

    ## Process papyrus
    # Get the Ids of all proteins that were processed via AlphaFold
    voc = VocSmiles()
    words = set()

    # Create output file

    out_file = os.path.join(output_dir, output_dataset_prefix + '.tsv')
    with open(out_file, 'w') as f:
        f.write('PID\tpchembl\ttokens\n')

    protein_id_idx = 9
    pchembl_val_idx = 21
    smiles_idx = 4

    # Processing line by line because of the dataset size
    with open(papyrus_file, 'r') as papyrus:
        print('Processing papyrus...')
        header = papyrus.readline()
        print([(i, attr) for i, attr in enumerate(header.split('\t'))])

        for i in tqdm(range(LENGTH_PAPYRUS)): # Process all papyrus entries
            try:
                entry = papyrus.readline()
                attributes = entry.split('\t')
                print()
                # If protein_id was processed by alphafold add an entry to new dataset
                protein_id = attributes[protein_id_idx]
                if protein_id in processed_proteins:
                    pchembl_val = attributes[pchembl_val_idx][:5]
                    smile = attributes[smiles_idx]

                    # Tokenize the smiles
                    tokenized = voc.split(smile)
                    if 10 < len(tokenized) <= 100: 
                        words.update(tokenized)
                        tokens = ' '.join(tokenized)

                        out = f'{protein_id}\t{pchembl_val}\t{tokens}\n'
                        open(out_file, 'a').write(out).close()
            except:
                print(entry)

    # Save vocabulary file
    if save_voc:
        print('Saving vocabulary...')
        with open(os.path.join(output_dir, 'voc_smiles.txt'), 'w') as voc:
            voc.write('\n'.join(sorted(words)))


def create_protein_subset_file(output_dir, write_out_file=False):
    """
    Parses the filtered papyrus files and creates a text file with all the unique
    protein accession numbers contained within it
    """

    target_id = 0

    directory = '/Users/amac/git/datasets'
    prefix = 'results_filtered_'
    names = ['ARs', 'ARs_high', 'Kin_human_high', 'CCRs_human_high', 'MonoamineR_human_high',
            'MonoamineR_human', 'SLC6_high']

    filenames = [prefix + suffix +  '.txt' for suffix in names]
    paths = [os.path.join(directory, f) for f in filenames]
    target = paths[target_id]
    out_file = os.path.join(output_dir, names[target_id] + '.txt')

    data = pd.read_csv(target, sep='\t')
    proteins = data['accession'].unique()

    # For writing out the output file that contains the filtered protein IDs
    if write_out_file:

        with open(out_file, 'w') as out:
            for entry in proteins:
                if '.' not in entry:
                    out.write(f'{entry}\n')



if __name__ == "__main__":

    af_output_dir = '/data/bernataviciusa/af_data/final'
    output_dir = '/data/bernataviciusa/af_data/data/kinases'
    papyrus_file = '/data/bernataviciusa/af_data/papyrus/papyrus_.tsv'
    output_dataset_prefix = 'dataset'
    protein_targets_file = '/data/bernataviciusa/af_data/papyrus/Kin_human_high.txt'


    build_dataset(af_output_dir, output_dir, papyrus_file, 
                  output_dataset_prefix, protein_targets_file=protein_targets_file, process=False)

    # create_subset('/Users/amac/git/datasets', write_out_file=True)

    # build_dataset()
    

