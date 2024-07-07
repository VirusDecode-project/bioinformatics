from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO
import subprocess

class SequenceAlignment:
    def __init__(self, files, muscle_exe="muscle"):
        self.files = files
        self.muscle_exe = muscle_exe
        self.combined_file = "combined.fasta"
        self.aligned_file = "aligned.fasta"
        self.aligned_dna_file = "aligned_dna.fasta"
        self.base_colors = {
            'A': '#ff9999',
            'T': '#99ff99',
            'G': '#9999ff',
            'C': '#ffff99',
            '-': '#cccccc'
        }
        self.sequences = []
        self.aligned_sequences = []
        self.protein_sequences = []
        self.mutations = []

    def read_sequences(self):
        self.sequences = [SeqIO.read(file, "fasta") for file in self.files]

    def combine_sequences(self):
        with open(self.combined_file, "w") as f:
            SeqIO.write(self.sequences, f, "fasta")

    def run_muscle(self):
        subprocess.run([self.muscle_exe, "-in", self.combined_file, "-out", self.aligned_file])

    def read_alignment(self):
        alignment = AlignIO.read(self.aligned_file, "fasta")
        self.aligned_sequences = [record.seq for record in alignment]

    def translate_sequences(self):
        self.protein_sequences = []
        for dna_seq in self.aligned_sequences:
            clean_dna_seq = str(dna_seq).replace('-', 'N')
            clean_dna_seq = Seq(clean_dna_seq)
            protein_seq = clean_dna_seq.translate()
            self.protein_sequences.append(protein_seq)
        for idx, protein_seq in enumerate(self.protein_sequences, start=1):
            print(f"Translated protein sequence {idx}: {protein_seq}")

    

    def run(self):
        self.read_sequences()
        self.combine_sequences()
        self.run_muscle()
        self.read_alignment()
        self.translate_sequences()


if __name__ == "__main__":
    files = ["NC_045512.fasta", "OX014251.fasta"]
    alignment = SequenceAlignment(files)
    alignment.run()
