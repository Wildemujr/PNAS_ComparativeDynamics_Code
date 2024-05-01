# Refactoring and Improving Imports

# Standard Libraries
import os
import sys
import csv
from pprint import pprint as pp
from ast import literal_eval
from itertools import product

# Third-party Libraries
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.SVDSuperimposer import SVDSuperimposer
from scipy.spatial.distance import cosine
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import warnings
from Bio import BiopythonWarning

# Suppress specific Biopython warnings
warnings.simplefilter('ignore', BiopythonWarning)


# Utility Functions

def read_fasta(file_path):
    """
    Reads sequences from a FASTA file.
    
    Parameters:
    - file_path (str): Path to the FASTA file.
    
    Returns:
    - list of str: List of sequences from the FASTA file.
    """
    with open(file_path, "r") as fasta_file:
        sequences = [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
    return sequences


def cosine_similarity(A, B):
    """
    Calculates the cosine similarity between two arrays.
    
    Parameters:
    - A (array-like): First array.
    - B (array-like): Second array.
    
    Returns:
    - float: Cosine similarity between A and B.
    """
    dot_product = np.dot(A, B)
    norm_A = np.linalg.norm(A)
    norm_B = np.linalg.norm(B)
    return dot_product / (norm_A * norm_B)


def frobenius_similarity(A, B):
    """
    Calculates the Frobenius norm of the difference between two matrices.
    
    Parameters:
    - A (matrix): First matrix.
    - B (matrix): Second matrix.
    
    Returns:
    - float: Frobenius norm of the difference between A and B.
    """
    diff = A - B
    return np.linalg.norm(diff, ord='fro')


def normalize_data(data):
    """
    Performs min-max normalization on the data.
    
    Parameters:
    - data (array-like): Data to be normalized.
    
    Returns:
    - array-like: Normalized data.
    """
    min_val = np.min(data)
    max_val = np.max(data)
    return (data - min_val) / (max_val - min_val)


def pearson_correlation(A, B):
    """
    Calculates the Pearson correlation coefficient between two arrays.
    
    Parameters:
    - A (array-like): First array.
    - B (array-like): Second array.
    
    Returns:
    - tuple: Pearson correlation coefficient and p-value.
    """
    return pearsonr(A, B)


def get_tiles(ANM_crossCorrelation_mat, sequence, lower_tile_length):
    """
    Extracts submatrices (tiles) and corresponding sequence tiles from a given matrix.
    
    Parameters:
    - ANM_crossCorrelation_mat: Numpy matrix from which tiles are extracted.
    - sequence: Corresponding sequence data.
    - lower_tile_length: Minimum length of tiles to extract.

    Returns:
    - List of dictionaries containing tiles and metadata.
    """
    N = len(ANM_crossCorrelation_mat)
    tiles = [
        {
            'tile': ANM_crossCorrelation_mat[start:end, start:end],
            'sequence': sequence[start:end],
            'length': L,
            'start_index': start + 1,
            'end_index': end,
            'center_index': np.median(np.arange(start, end)) + 1
        }
        for L in range(N, lower_tile_length-1, -1)
        for start in range(N - L + 1)
        for end in [start + L]
    ]
    return tiles


def get_residue_range(pdb_file, chain_id):
    """
    Returns the range of residue numbers for a given chain in a PDB file.
    
    Parameters:
    - pdb_file: Path to the PDB file.
    - chain_id: ID of the chain for which the residue range is to be found.

    Returns:
    - Tuple of (minimum residue number, maximum residue number).
    """
    # Initialize a PDBParser object
    parser = PDBParser()

    # Read the structure
    structure = parser.get_structure('pdb_structure', pdb_file)

    # Get the specific chain
    chain = structure[0][chain_id]

    # Extract residue numbers and determine min and max
    residue_numbers = [residue.get_id()[1] for residue in chain]
    return min(residue_numbers), max(residue_numbers)


def extract_coordinates(pdb_file, chain_id, start, end):
    """
    Extracts the alpha carbon (CA) coordinates from a PDB file for a given residue range.

    Parameters:
    - pdb_file (str): Path to the PDB file.
    - chain_id (str): Chain ID to extract coordinates from.
    - start (int): Start residue number.
    - end (int): End residue number.

    Returns:
    - List of tuples: A list of 3D coordinates.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    chain = structure[0][chain_id]

    coords = []

    for residue in chain:
        if start <= residue.get_id()[1] <= end:
            ca_atom = residue['CA']
            coords.append(ca_atom.get_coord())

    return coords


def group_tiles_by_length(tiles):
    """
    Groups tiles by their length.
    
    Parameters:
    - tiles: List of tiles.

    Returns:
    - Dictionary with keys as tile lengths and values as lists of tiles of that length.
    """
    tiles_by_length = {}
    for index, tile in tiles.iterrows():
        tiles_by_length.setdefault(tile['length'], []).append(tile)
    
    return tiles_by_length


def write_similarity_to_csv(writer, tile_1, tile_2, length, similarity, frobenius_sim, within_protein=False):
    """
    Writes the similarity data to a CSV file. Can handle both within-protein and between-proteins comparisons.
    
    Parameters:
    - writer: CSV writer object.
    - tile_1: Tile from the first protein or first tile from the same protein.
    - tile_2: Tile from the second protein or second tile from the same protein.
    - length: Length of the tile.
    - similarity: Calculated cosine similarity.
    - frobenius_sim: Calculated Frobenius similarity.
    - within_protein: Boolean flag to indicate if the comparison is within the same protein.
    """
    if within_protein:
        writer.writerow([
            # tile_1['sequence'], tile_1['start_index'], tile_1['end_index'], tile_1['center_index'],
            # tile_2['sequence'], tile_2['start_index'], tile_2['end_index'], tile_2['center_index'],
            # length, similarity, frobenius_sim
            tile_1['sequence'], tile_2['sequence'], 
            tile_1['start_index'], tile_1['end_index'], tile_1['center_index'], 
            tile_2['start_index'], tile_2['end_index'], tile_2['center_index'], 
            length, similarity, frobenius_sim
        ])
    else:
        writer.writerow([
            tile_1['sequence'], tile_2['sequence'], 
            tile_1['start_index'], tile_1['end_index'], tile_1['center_index'], 
            tile_2['start_index'], tile_2['end_index'], tile_2['center_index'], 
            length, similarity, frobenius_sim
        ])


def output_superimposed_pdb(P, Q_rotated, iteration_number):
    """Output a PDB file with two superimposed structures."""
    
    # Convert your numpy arrays P and Q_rotated to Bio.PDB structures if they aren't already
    # This step will vary based on how your data is structured
    
    io = PDBIO()
    io.set_structure(P)  # Assuming P is a Bio.PDB structure
    io.save(f"superimposed_iteration_{iteration_number}_P.pdb")
    
    io.set_structure(Q_rotated)  # Assuming Q_rotated is a Bio.PDB structure
    io.save(f"superimposed_iteration_{iteration_number}_Q.pdb")



def kabsch_rmsd(P, Q, threshold=2.0, max_iterations=15, convergence_ratio=0.7):
    """
    Calculate the RMSD between two sets of 3D coordinates using the Kabsch algorithm 
    with iterative outlier rejection.

    Parameters:
    - P (array-like): First set of 3D coordinates.
    - Q (array-like): Second set of 3D coordinates.
    - threshold (float): Distance threshold for outlier rejection.
    - max_iterations (int): Maximum number of iterations for outlier rejection.
    - convergence_ratio (float): The ratio below which if the remaining points fall, the iteration stops.

    Returns:
    - float: RMSD value.
    """

    original_length = len(P)
    
    for iteration in range(max_iterations):
        # Translate points to their centroids
        centroid_P = np.mean(P, axis=0)
        centroid_Q = np.mean(Q, axis=0)
        P -= centroid_P
        Q -= centroid_Q

        # Compute the cross-covariance matrix
        H = np.dot(P.T, Q)

        # Compute the rotation matrix using SVD
        U, _, Vt = np.linalg.svd(H)
        R = np.dot(Vt.T, U.T)

        # Handle reflections
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = np.dot(Vt.T, U.T)

        # Rotate P
        P_rotated = np.dot(P, R)

        # Identify outliers
        distances = np.sqrt(np.sum((P_rotated - Q) ** 2, axis=1))
        outliers = distances > threshold

        # Check for convergence condition
        if np.sum(outliers) / original_length > 1 - convergence_ratio:
            break

        # Remove outliers
        P = P[~outliers]
        Q = Q[~outliers]

    # Compute the final RMSD
    rmsd = np.sqrt(np.mean((P_rotated - Q) ** 2))
    return rmsd




def calculate_cosine_similarities_between_proteins(tiles_protein_1, tiles_protein_2, modeNum, pdb_file_1, pdb_file_2, chain_id_1, chain_id_2):
    """
    Calculates the cosine similarities between tiles of two proteins and saves the results to a CSV file.
    
    Parameters:
    - tiles_protein_1: List of tiles from the first protein.
    - tiles_protein_2: List of tiles from the second protein.
    - modeNum: Mode number (used for naming the output CSV file).
    - pdb_file_1: Path to the PDB file of the first protein.
    - pdb_file_2: Path to the PDB file of the second protein.
    - chain_id_1: Chain ID for the first protein.
    - chain_id_2: Chain ID for the second protein.
    """
    # Group tiles by their length
    tiles_by_length_1 = group_tiles_by_length(tiles_protein_1)
    tiles_by_length_2 = group_tiles_by_length(tiles_protein_2)

    sup = SVDSuperimposer()

    with open(f'cosine_similarities_Mode{modeNum}.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["sequence_1", "sequence_2", "start_index_1", "end_index_1", "center_index_1", 
                         "start_index_2", "end_index_2", "center_index_2", "tile_length", 
                         "cosine_similarity", "frobenius_similarity"])


        # Iterate over lengths present in both proteins
        for length in set(tiles_by_length_1.keys()) & set(tiles_by_length_2.keys()):
            for tile_2 in tqdm(tiles_by_length_2[length], desc=f"Processing length {length}"):
                for tile_1 in tiles_by_length_1[length]:
                    flat_tile_1 = tile_1['tile'].flatten()
                    flat_tile_2 = tile_2['tile'].flatten()
                    similarity = cosine_similarity(flat_tile_1, flat_tile_2)
                    frobenius_sim = frobenius_similarity(tile_1['tile'], tile_2['tile'])
                    
                    # Write the results to the CSV file
                    write_similarity_to_csv(writer, tile_1, tile_2, length, similarity, frobenius_sim)

    print("Finished comparing all matching tile lengths.")



def calculate_cosine_similarities_within_protein(tiles_protein, modeNum):
    """
    Calculates the cosine similarities between tiles within the same protein and saves the results to a CSV file.
    
    Parameters:
    - tiles_protein: List of tiles from the protein.
    - modeNum: Mode number (used for naming the output CSV file).
    """
    # Group tiles by their length
    tiles_by_length = group_tiles_by_length(tiles_protein)

    with open(f'cosine_similarities_Mode{modeNum}_within_protein.csv', 'w', newline='') as file:
        writer = csv.writer(file)

        writer.writerow(["sequence_1", "sequence_2", "start_index_1", "end_index_1", "center_index_1", 
                    "start_index_2", "end_index_2", "center_index_2", "tile_length", 
                    "cosine_similarity", "frobenius_similarity"])
        
        # writer.writerow(["sequence_1", "start_index_1", "end_index_1", "center_index_1", 
        #                  "sequence_2", "start_index_2", "end_index_2", "center_index_2", 
        #                  "tile_length", "cosine_similarity", "frobenius_similarity"])

        # Iterate over all lengths
        for length, tiles in tiles_by_length.items():
            for tile_1 in tqdm(tiles, desc=f"Processing length {length}"):
                for tile_2 in tiles:
                    flat_tile_1 = tile_1['tile'].flatten()
                    flat_tile_2 = tile_2['tile'].flatten()
                    similarity = cosine_similarity(flat_tile_1, flat_tile_2)
                    frobenius_sim = frobenius_similarity(tile_1['tile'], tile_2['tile'])
                    write_similarity_to_csv(writer, tile_1, tile_2, length, similarity, frobenius_sim, within_protein=True)

    print("Finished comparing all matching tile lengths within the protein.")



# Here's a refactored version that breaks down the main function:

def read_and_process_cross_correlation_matrix(csv_dir, pdb_id, anm_mode_number, first_res, last_res):
    """
    Reads and processes the cross-correlation matrix.
    
    Parameters:
    - csv_dir: Directory containing the CSV file.
    - pdb_id: ID of the protein.
    - anm_mode_number: Mode number for ANM.
    - first_res: First residue index.
    - last_res: Last residue index.

    Returns:
    - Processed cross-correlation matrix.
    """
    csv_file = f"{pdb_id}_CrossCorr_Mode{anm_mode_number}.csv"
    csv_path = os.path.join(csv_dir, csv_file)
    
    ANM_crossCorrelation_mat = pd.read_csv(csv_path, header=None).values
    return ANM_crossCorrelation_mat[first_res:last_res + 1, first_res:last_res + 1]




def process_protein_tiles(ANM_crossCorrelation_mat, sequence, pdb_id):
    """
    Processes the protein tiles. If the pickle file exists, it reads it. 
    Otherwise, it computes the tiles and saves them as a pickle file.
    
    Parameters:
    - ANM_crossCorrelation_mat: The cross-correlation matrix.
    - sequence: Protein sequence.
    - pdb_id: ID of the protein.

    Returns:
    - DataFrame of protein tiles.
    """
    pickle_file = f'{pdb_id}_tiles.pkl'
    
    if not os.path.exists(pickle_file):
        tiles_protein = get_tiles(ANM_crossCorrelation_mat, sequence, 6)
        prot_tile_df = pd.DataFrame.from_dict(tiles_protein).sort_values(by=['length', 'center_index']).reset_index(drop=True)
    else:
        prot_tile_df = pd.read_pickle(pickle_file)
    
    return prot_tile_df



def process_output_csv(anm_mode_number, protein_length, within_protein=False):
    """
    Processes the output CSV file. Calculates and adds the normalized frobenius similarity.
    Filters the data based on tile length and cosine similarity conditions.

    Parameters:
    - anm_mode_number: Mode number for ANM.
    - protein_length: Length of the protein sequence.

    Returns:
    - DataFrame containing the filtered average cosine similarities.
    """
    
    if within_protein:
        csv_file = f'cosine_similarities_Mode{anm_mode_number}_within_protein.csv'
    else:
        csv_file = f'cosine_similarities_Mode{anm_mode_number}.csv'
    
    tile_df = pd.read_csv(csv_file)
    tile_df['frobenius_similarity_normed'] = normalize_data(tile_df['frobenius_similarity'])
    tile_df.to_csv(csv_file, index=False)

    threshold = 1e-9
    condition = (tile_df['tile_length'] >= 6) & (tile_df['tile_length'] <= protein_length) & (np.abs(tile_df['cosine_similarity'] - 1) > threshold)
    filtered_tile_df = tile_df[condition]
    
    grouped = filtered_tile_df.groupby(['center_index_1', 'tile_length'])
    average_cosine_similarity = grouped['cosine_similarity'].mean()
    counts = grouped.size()

    result = pd.DataFrame({'average_cosine_similarity': average_cosine_similarity, 'counts': counts}).reset_index()
    result.rename(columns={
        'center_index_1': 'Tile_Center_Index',
        'tile_length': 'Tile_Length',
        'average_cosine_similarity': 'Average_Cosine_Similarity',
        'counts': 'Count'
    }, inplace=True)
    
    return result








# The main function is then simplified as follows:

def main():
    args = sys.argv
    csv_dir = os.path.join("../", "ANM_modes")
    print(args)
    
    if len(args) == 4:
        pdb_id, anm_mode_number, chain_id = args[1:4]

        pdb_file = f"{pdb_id}.chain{chain_id}.nowa.pdb"
        first_res, last_res = get_residue_range(pdb_file, chain_id)
        
        ANM_crossCorrelation_mat = read_and_process_cross_correlation_matrix(csv_dir, pdb_id, anm_mode_number, first_res, last_res)
        concat_sequence = read_fasta(os.path.join(os.getcwd(), f"{pdb_id}.nowa.fasta"))
        prot_tile_df = process_protein_tiles(ANM_crossCorrelation_mat, concat_sequence[0], pdb_id)
        
        protein_length = len(concat_sequence[0])
        calculate_cosine_similarities_within_protein(prot_tile_df, anm_mode_number)
        
        result = process_output_csv(anm_mode_number, protein_length, within_protein = True)
        result.to_csv(f'filtered_average_cosine_similarities_Mode{anm_mode_number}_within_protein.csv', index=False)
    
    elif len(args) == 5 or len(args) == 6:
        pdb_id_one, pdb_id_two, anm_mode_number, chain_id_one, chain_id_two = args[1:6]

        # Process for the first protein
        pdb_file_one = f"{pdb_id_one}.chain{chain_id_one}.nowa.pdb"
        first_res_one, last_res_one = get_residue_range(pdb_file_one, chain_id_one)
        ANM_crossCorrelation_mat_1 = read_and_process_cross_correlation_matrix(csv_dir, pdb_id_one, anm_mode_number, first_res_one, last_res_one)
        concat_sequence_A = read_fasta(os.path.join(os.getcwd(), f"{pdb_id_one}.chain{chain_id_one}.nowa.fasta"))
        prot_tile_df_one = process_protein_tiles(ANM_crossCorrelation_mat_1, concat_sequence_A[0], pdb_id_one)
        
        # Process for the second protein
        pdb_file_two = f"{pdb_id_two}.nowa.pdb"
        first_res_two, last_res_two = get_residue_range(pdb_file_two, chain_id_two)
        ANM_crossCorrelation_mat_2 = read_and_process_cross_correlation_matrix(csv_dir, pdb_id_two, anm_mode_number, first_res_two, last_res_two)
        concat_sequence_B = read_fasta(os.path.join(os.getcwd(), f"{pdb_id_two}.nowa.fasta"))
        prot_tile_df_two = process_protein_tiles(ANM_crossCorrelation_mat_2, concat_sequence_B[0], pdb_id_two)
        
        # Calculate the cosine similarities between the two proteins
        calculate_cosine_similarities_between_proteins(prot_tile_df_one, prot_tile_df_two, anm_mode_number, pdb_file_one, pdb_file_two, chain_id_one, chain_id_two)
        min_protein_length = min(len(concat_sequence_A[0]), len(concat_sequence_B[0]))
        result = process_output_csv(anm_mode_number, min_protein_length, within_protein = False)
        result.to_csv(f'filtered_average_cosine_similarities_Mode{anm_mode_number}_between_chain_{chain_id_one}_and_{chain_id_two}.csv', index=False)


if __name__ == "__main__":
    main()