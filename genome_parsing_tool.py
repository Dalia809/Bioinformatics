import argparse
from Bio.Seq import Seq

def read_fasta(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

def find_open_reading_frames(sequence, strand='forward'):
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    for frame in range(3):
        start = None
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon == 'ATG' and start is None:
                start = i
            elif codon in stop_codons and start is not None:
                orfs.append((sequence[start:i+3], start, i+3, strand, frame))
                start = None
    return orfs

def main():
    parser = argparse.ArgumentParser(description="Extract regions between start and stop codons from a FASTA file, including reverse complement.")
    parser.add_argument('input_file', help="Path to the input FASTA file containing the genome sequence")
    args = parser.parse_args()
    sequence = read_fasta(args.input_file)
    orfs_forward = find_open_reading_frames(sequence, strand='forward')
    # Used Chat gpt to use the biopython library to implement the reverse complement logic using the following prompt:
    # Use biopython to extend this tool to include the reverse complement and search in all six possible reading frames for genes. {Pasted code From first commit here}
    reverse_complement_seq = str(Seq(sequence).reverse_complement())
    orfs_reverse = find_open_reading_frames(reverse_complement_seq, strand='reverse')
    orfs = orfs_forward + orfs_reverse
    for i, (orf, start, end, strand, frame) in enumerate(orfs, 1):
        print(f"ORF {i}: {orf}, Start: {start}, End: {end}, Strand: {strand}, Frame: {frame + 1}")

if __name__ == '__main__':
    main()
