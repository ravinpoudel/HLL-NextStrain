import os
import pickle
from Bio import SeqIO
from datasketch import HyperLogLog

def extract_kmers(sequence, k):
    """
    Extract k-mers from a given sequence.

    Args:
        sequence (str): The input nucleotide sequence.
        k (int): The length of k-mers.

    Returns:
        list: A list of k-mers.
    """
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]

def create_hll_from_sequence(sequence, k=31):
    """
    Create a HyperLogLog (HLL) data structure from a given sequence.

    Args:
        sequence (str): The input nucleotide sequence.
        k (int, optional): The length of k-mers. Defaults to 31.

    Returns:
        HyperLogLog: The HLL representation of the sequence.
    """
    hll = HyperLogLog()
    kmers = extract_kmers(sequence, k)
    for kmer in kmers:
        hll.update(kmer.encode('utf8'))
    return hll

def load_existing_hlls(hll_dir):
    """
    Load existing HLLs from files in a specified directory.

    Args:
        hll_dir (str): The directory containing HLL files.

    Returns:
        dict: A dictionary with strain names as keys and HLLs as values.
    """
    hlls = {}
    for filename in os.listdir(hll_dir):
        if filename.endswith(".hll"):
            strain_name = filename.replace(".hll", "")
            with open(os.path.join(hll_dir, filename), "rb") as f:
                hlls[strain_name] = pickle.load(f)
    return hlls

def save_hll(hll, filename):
    """
    Save an HLL to a file.

    Args:
        hll (HyperLogLog): The HLL to save.
        filename (str): The file to save the HLL to.
    """
    with open(filename, "wb") as f:
        pickle.dump(hll, f)

def compare_with_existing_strains(new_hll, existing_hlls):
    """
    Compare a new strain's HLL with existing strains' HLLs.

    Args:
        new_hll (HyperLogLog): The HLL of the new strain.
        existing_hlls (dict): A dictionary with strain names as keys and HLLs as values.

    Returns:
        dict: A dictionary with strain names as keys and Jaccard similarities as values.
    """
    similarities = {}
    for strain_name, hll in existing_hlls.items():
        jaccard_similarity = new_hll.jaccard(hll)
        similarities[strain_name] = jaccard_similarity
    return similarities

def main():
    """
    Main function to create HLLs for existing strains, create an HLL for a new strain,
    and compare the new strain with existing strains.
    """
    existing_strains_dir = "./existing_strains"  # Directory containing existing strain sequences in FASTA format
    hll_dir = "./hll_data"  # Directory to save/load HLLs
    new_strain_file = "./new_strain.fasta"  # New strain sequence in FASTA format

    # Step 1: Create and save HLLs for existing strains (run once)
    os.makedirs(hll_dir, exist_ok=True)
    for filename in os.listdir(existing_strains_dir):
        if filename.endswith(".fasta"):
            strain_name = filename.replace(".fasta", "")
            file_path = os.path.join(existing_strains_dir, filename)
            with open(file_path, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    hll = create_hll_from_sequence(str(record.seq))
                    save_hll(hll, os.path.join(hll_dir, f"{strain_name}.hll"))

    # Step 2: Load HLLs of existing strains
    existing_hlls = load_existing_hlls(hll_dir)

    # Step 3: Create HLL for the new strain
    with open(new_strain_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            new_hll = create_hll_from_sequence(str(record.seq))

    # Step 4: Compare the new strain with existing strains
    similarities = compare_with_existing_strains(new_hll, existing_hlls)
    for strain_name, similarity in similarities.items():
        print(f"Jaccard similarity with {strain_name}: {similarity}")

if __name__ == "__main__":
    main()

