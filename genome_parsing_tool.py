import argparse
from Bio.Seq import Seq
import os

def read_fasta(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

def find_orfs(sequence):
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    
    for frame in range(3):
        i = frame
        while i < len(sequence) - 2:
            codon = sequence[i:i+3]
            if codon == 'ATG': 
                start = i
                for j in range(i + 3, len(sequence) - 2, 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in stop_codons:
                        orfs.append(sequence[start:j+3])
                        break
            i += 3
    return orfs

def find_all_orfs(sequence):
    orfs = []
    orfs.extend(find_orfs(sequence)) 
    
    reverse_complement_seq = str(Seq(sequence).reverse_complement())
    orfs.extend(find_orfs(reverse_complement_seq))  
    
    return orfs

def translate_orfs_to_proteins(orfs, min_length):
    proteins = set()
    for orf in orfs:
        protein = str(Seq(orf).translate(to_stop=True))
        if len(protein) >= min_length:
            proteins.add((protein, len(protein)))
    return proteins

def process_file(file_path, min_length):
    sequence = read_fasta(file_path)
    orfs = find_all_orfs(sequence)
    proteins = translate_orfs_to_proteins(orfs, min_length)
    
    file_dir = os.path.dirname(file_path)
    file_name = os.path.basename(file_path)
    
    output_file = os.path.join(file_dir, f"{file_name}_output.txt")
    
    with open(output_file, 'w') as f:
        f.write(f"Proteins for file {file_path}:\n")
        for protein, length in proteins:
            f.write(f"{protein} {length}\n")
        f.write("\n")

def main():
    parser = argparse.ArgumentParser(description="Extract ORFs, filter by length, and translate them into protein strings.")
    parser.add_argument('input_files', nargs='+', help="Paths to the input FASTA files containing the DNA sequences")
    parser.add_argument('--min_length', type=int, default=100, help="Minimum length of ORFs to be considered (in codons, default: 100)")
    
    args = parser.parse_args()
    
    for file_path in args.input_files:
        if "GCA" in file_path:
            process_file(file_path, args.min_length)
main()
