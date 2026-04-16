import numpy as np
import pandas as pd
from cyvcf2 import VCF
from scipy.sparse import coo_matrix, save_npz, load_npz
from collections import defaultdict
from tqdm import tqdm
import sys
import os
import pickle

DATA_FILEPATHS = {
    'manifest': '../illumina_manifests/GSA-manifest-build37.csv',
    'target_vcf': '../chromosomes/HG00096-chr13.vcf.gz', # any chr13 vcf will do
    'source_snp': 'rs9315973',
    'source_snp_id': 'rs9315973',
    'source_snp_chromosome': '13',
    'max_ld_association_depth': 15
}

def get_or_create_node_id(snp_str):
    """
    Return a serialized SNP ID from the map and a serialized ID for the next
    unseen SNP to be inserted into the map.
    """
    global next_id, SnpIdMap
    if snp_str not in SnpIdMap:
        SnpIdMap[snp_str] = next_id
        next_id += 1
    return SnpIdMap[snp_str]

def read_snp_pairs(filename):
    """
    Generator yielding (snp_a, snp_b) tuples from each valid line of the file.
    """
    with open(filename, "r") as f:
        for line in f:
            # Skip header or empty lines
            if "CHR_A" in line or not line.strip():
                continue
            # fields: CHR_A, BP_A, SNP_A, MAF_A, CHR_B, BP_B, SNP_B, MAF_B, R2 
            parts = line.split()
            yield parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6], parts[-1]

def get_neighbors(adjacency_mat, node_id, threshold=0.5):
    """
    Given an adjacency matrix and a serialized node ID, return the column 
    indices of neighbors with weights above a specified threshold.
    """
    neighbors = adjacency_mat[node_id].nonzero()[1]
    weights = adjacency_mat[node_id].data
    curr_col_weight_pairs = [(col, w) for col, w in zip(neighbors, weights) if w > threshold]
    neighbors_idx = [pair[0] for pair  in curr_col_weight_pairs]
    return neighbors_idx

if __name__ == '__main__':
    
    df = pd.read_csv(DATA_FILEPATHS['manifest'])
    df = df[df['Chr'].notna()]
    df = df[['IlmnID', 'Name', 'SNP', 'Chr', 'MapInfo']]
    df['MapInfo'] = df['MapInfo'].astype(int)
    df['Chr'] = df['Chr'].astype(str)
    genotyped_chromosome = df[df['Chr'] == DATA_FILEPATHS['source_snp_chromosome']]
    
    vcf_filepath = DATA_FILEPATHS['target_vcf']
    vcf = VCF(vcf_filepath, gts012=True)

    meta_info = {'CHROM':[], 'POS':[], 'REF':[], 'ALT':[], 'ID': [], 'AlleleFreq': []}
    for v in tqdm(vcf):
        meta_info['CHROM'].append(v.CHROM)
        meta_info['POS'].append(v.POS)
        meta_info['REF'].append(v.REF)
        meta_info['ALT'].append(v.ALT[0])
        meta_info['ID'].append(v.ID)
        meta_info['AlleleFreq'].append(v.INFO['AF'])
        
    meta = pd.DataFrame(meta_info)
    
    global next_id, SnpIdMap
    SnpIdMap = {}
    next_id = 0
    counter = 0
    neighbor_count = defaultdict(int)
    if len(sys.argv) != 2:
        print("Usage: python3 build_LD_network.py <file_to_process> (PLINK output with pairwise LD associations)")
        sys.exit(1)
    file_to_process = sys.argv[1]
    
    if not os.path.exists(f"{file_to_process.split('-')[1]}_network_weighted_adjacency.npz"):
        # establish the SNP to serialized ID mapping
        for chr_a, _, snp_a, _, chr_b, _, snp_b, r2 in read_snp_pairs(file_to_process):
            id_a = get_or_create_node_id(snp_a)
            id_b = get_or_create_node_id(snp_b)
            counter += 1
            if counter % 10000000 == 0:
                print(f'{counter} pairwise relationships processed')
            if chr_a == 'X' or chr_a == 'x':
                print(snp_a, snp_b, r2)
                break
        # build adjacency matrix
        SnpIdMap_reversed = {v: k for k, v in SnpIdMap.items()}
        n_rows = len(SnpIdMap)
        n_cols = n_rows
        row_idx = []
        col_idx = []
        edge_weights = []
        counter = 0
        for _, _, snp_a, _, _, _, snp_b, r2 in read_snp_pairs(file_to_process):
            row_num = SnpIdMap[snp_a]
            col_num = SnpIdMap[snp_b]
            
            row_idx.append(row_num)
            col_idx.append(col_num)
            edge_weights.append(float(r2))
            
            row_idx.append(col_num)
            col_idx.append(row_num)
            edge_weights.append(float(r2))
            
            counter += 1
            if counter % 10000000 == 0:
                print(f'{counter} pairwise relationships inserted to adjacency matrix')
        network_coo_matrix = coo_matrix(
            (edge_weights, (row_idx, col_idx)), 
            shape=(n_rows, n_cols), 
            dtype=np.float64)
        network_weighted_adjacency = network_coo_matrix.tocsr()
        save_npz(f"{file_to_process.split('-')[1]}_network_weighted_adjacency.npz", network_weighted_adjacency)
        with open(f"{file_to_process.split('-')[1]}_snp_id_serialized_mapping.pkl", 'wb') as f:
            pickle.dump(SnpIdMap, f)
    else:
        network_weighted_adjacency = load_npz(f"{file_to_process.split('-')[1]}_network_weighted_adjacency.npz")
        with open(f"{file_to_process.split('-')[1]}_snp_id_serialized_mapping.pkl", 'rb') as f:
            SnpIdMap = pickle.load(f)
        SnpIdMap_reversed = {v: k for k, v in SnpIdMap.items()}
            
    # target SNP to find the n-th degree neighbors for
    ##############################################################
    source_snp = SnpIdMap[DATA_FILEPATHS['source_snp_id']]
    ##############################################################
    min_ld_threshold = float(file_to_process.split('-')[-1][:-4])
    neighbor_degree_mapping = {}
    # adjust the min LD parameter in DATA_FILEPATHS if needed, all edge 
    # associations with r2 > specified are treated as (connected) neighbors
    max_degree = DATA_FILEPATHS['max_ld_association_depth']
    all_k_degree_neighbors = set([source_snp])
    neighbor_degree_mapping[0] = set([DATA_FILEPATHS['source_snp_id']])
    directly_genotyped = set(meta[meta['POS'].isin(genotyped_chromosome['MapInfo'])]['ID'].to_list())
    
    for iteration in range(1, max_degree + 1):
        all_k_degree_neighbors_deepcopy = set(all_k_degree_neighbors)
        for node in all_k_degree_neighbors:
            curr_node_neighbors = get_neighbors(
                network_weighted_adjacency, node, threshold=min_ld_threshold
            )
            all_k_degree_neighbors_deepcopy.update(curr_node_neighbors)
        
        curr_degree_snps = all_k_degree_neighbors_deepcopy.difference(all_k_degree_neighbors)
        curr_degree_snps_ids_unserialized = [SnpIdMap_reversed[int(i)] for i in list(curr_degree_snps)]
        target_neighbor_chromosome_meta = meta[meta['ID'].isin(pd.Series(curr_degree_snps_ids_unserialized))]
        target_neighbor_genotyped = target_neighbor_chromosome_meta[
            target_neighbor_chromosome_meta['POS'].isin(genotyped_chromosome['MapInfo'])]

        neighbor_degree_mapping[iteration] = set(target_neighbor_genotyped['ID'].to_list())
        all_k_degree_neighbors = set(all_k_degree_neighbors_deepcopy)

        k_degree_neighbor_ids_unserialized = []
        for idx in all_k_degree_neighbors:
            k_degree_neighbor_ids_unserialized.append(SnpIdMap_reversed[int(idx)])
        # make a directory for the snp first if there isn't one
        with open(f"./{DATA_FILEPATHS['source_snp']}-neighbors-by-degree/{iteration}-degree-neighbors-r2-{min_ld_threshold}.txt", 'w') as f:
            for snp in k_degree_neighbor_ids_unserialized:
                if snp in directly_genotyped:
                    f.write(f"{snp}\n")
        print(f'{len(k_degree_neighbor_ids_unserialized)} neighbors in {iteration}-th degree neighborhood with r2 > {min_ld_threshold}')
    
    total = 0
    for i in neighbor_degree_mapping.items():
        total += len(i[1])
        print(f'degree {i[0]}: {len(i[1])}, total: {total}')

    with open(f"{DATA_FILEPATHS['source_snp']}_neighbor_degree_mapping.pkl", 'wb') as f:
        pickle.dump(neighbor_degree_mapping, f)
        