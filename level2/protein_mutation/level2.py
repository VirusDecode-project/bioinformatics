from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO, SeqIO
import subprocess

class SequenceAlignment:
    def __init__(self, files, reference_file, muscle_exe="muscle"):
        self.reference_file = reference_file
        self.files = files
        self.files_len = len(files)
        self.muscle_exe = muscle_exe
        self.combined_file = "combined.fasta"
        self.aligned_file = "aligned.fasta"
        self.protein_file = "proteins.fasta"
        self.aligned_protein_file = "aligned_proteins.fasta"
        self.reference_sequence = None
        self.sequences = []
        self.aligned_sequences = []
        self.protein_sequences = {}
        self.mutations = []
        self.reference_index = None
        self.reference_id = None
        self.variants_id =[]

    def read_sequences(self):
        self.reference_sequence = SeqIO.read(self.reference_file, "fasta")
        self.sequences = [SeqIO.read(file, "fasta") for file in self.files]
        self.reference_id = self.reference_sequence.id
        self.variants_id = [record.id for record in self.sequences]


    def combine_sequences(self):
        with open(self.combined_file, "w") as f:
            SeqIO.write([self.reference_sequence] + self.sequences, f, "fasta")

    def run_muscle_dna(self):
        result = subprocess.run([self.muscle_exe, "-in", self.combined_file, "-out", self.aligned_file])
        if result.returncode != 0:
            print("Error running muscle:")
            print(result.stderr)
        else:
            print("Muscle ran successfully:")
            print(result.stdout)

    def read_alignment(self):
        alignment = AlignIO.read(self.aligned_file, "fasta")
        self.aligned_sequences = [record for record in alignment]
        
        print(self.reference_sequence.id)
        for i, record in enumerate(self.aligned_sequences):
            if record.id == self.reference_sequence.id:
                self.reference_index = i
                break

    def translate_sequences(self):
        for record in self.aligned_sequences:
            clean_dna_seq = str(record.seq).replace('-', 'N')
            clean_dna_seq = Seq(clean_dna_seq)
            protein_seq = clean_dna_seq.translate()
            self.protein_sequences[record.id] = SeqRecord(protein_seq, id=record.id)

    def write_protein_sequences(self):
        with open(self.protein_file, "w") as f:
            SeqIO.write(self.protein_sequences.values(), f, "fasta")

    def find_mutation(self):
        for variant_id in self.variants_id:
            if self.reference_index is not None:
                mutation = [
                    i for i, (a, b) in enumerate(zip(self.protein_sequences[self.reference_id].seq, self.protein_sequences[variant_id].seq))
                    if a != b and b != 'X'
                ]
                self.mutations.append((variant_id, mutation))
                print(f"Total mutations for variant {variant_id}: {len(mutation)}/{len(self.protein_sequences[self.reference_id].seq)}")
            else:
                print("Reference sequence not found in alignment.")

    def print_mutation(self):
        for variant_id, mutation in self.mutations:
            print(f"Mutations for variant {variant_id}:")
            for mut in mutation:
                print(f"Position {mut + 1}: {self.protein_sequences[self.reference_id].seq[mut]} -> {self.protein_sequences[variant_id].seq[mut]}")

    def run(self):
        self.read_sequences()
        self.combine_sequences()
        self.run_muscle_dna()
        self.read_alignment()
        self.translate_sequences()
        self.write_protein_sequences()
        self.find_mutation()
        self.print_mutation()

if __name__ == "__main__":
    reference_file = "reference.fasta"
    files = ["variant2.fasta", "variant3.fasta"]
    alignment = SequenceAlignment(files, reference_file)
    alignment.run()
