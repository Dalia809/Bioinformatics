import argparse
from Bio.Seq import Seq

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

# Used chat gpt to write the translation function from orfs to protiens and the prompt was:
# Write me a function that accepts the orfs and translates them to protiens and returns them
def translate_orfs_to_proteins(orfs):
    proteins = set()
    for orf in orfs:
        protein = str(Seq(orf).translate(to_stop=True))
        if protein:
            proteins.add(protein)
    return proteins

def main():
    parser = argparse.ArgumentParser(description="Extract ORFs and translate them into protein strings.")
    parser.add_argument('input_file', help="Path to the input FASTA file containing the DNA sequence")
    args = parser.parse_args()
    
    sequence = read_fasta(args.input_file)
    
    orfs = find_all_orfs(sequence)    
    proteins = translate_orfs_to_proteins(orfs)
    for protein in proteins:
        print(protein)

if __name__ == '__main__':
    main()
