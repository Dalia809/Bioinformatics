import argparse

def read_fasta(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

def find_open_reading_frames(sequence):
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    for frame in range(3):
        start = None
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon == 'ATG' and start is None:
                start = i
            elif codon in stop_codons and start is not None:
                orfs.append(sequence[start:i+3])
                start = None
    return orfs

def main():
    parser = argparse.ArgumentParser(description="Extract regions between start and stop codons from a FASTA file.")
    parser.add_argument('input_file', help="Path to the input FASTA file containing the genome sequence")
    
    args = parser.parse_args()
    sequence = read_fasta(args.input_file)
    orfs = find_open_reading_frames(sequence)
    for i, orf in enumerate(orfs, 1):
        print(f"ORF {i}: {orf}")

main()