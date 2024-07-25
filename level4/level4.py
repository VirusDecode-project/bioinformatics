from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO, SeqIO, Entrez
from urllib.error import HTTPError
import subprocess



class SequenceAlignment:
    def __init__(self, files, reference_id, muscle_exe="muscle"):
        Entrez.email = "your_email@example.com"
        self.files = files
        self.muscle_exe = muscle_exe
        self.combined_file = "combined.fasta"
        self.aligned_file = "aligned.fasta"
        self.protein_file = "proteins.fasta"
        self.reference_sequence = None
        self.variant_sequences = []
        self.aligned_sequences = []
        self.protein_sequences = {}
        self.mutations = []
        self.reference_index = None
        self.reference_id = reference_id
        self.CDS_dict = {}

    def read_sequences(self):
        try:
            # Get reference sequence
            handle = Entrez.efetch(db="nucleotide", id=reference_id, rettype="gb", retmode="text")
            self.reference_sequence = SeqIO.read(handle, "genbank")
            handle.close()

            # Get variant sequences
            self.variant_sequences = [SeqIO.read(file, "fasta") for file in self.files]
            
            # Combine reference and variant sequences
            with open(self.combined_file, "w") as f:
                SeqIO.write([self.reference_sequence] + self.variant_sequences, f, "fasta")

            self.variants_id = [record.id for record in self.variant_sequences]
        except HTTPError as e:
            print(f"HTTPError: {e.code} - {e.reason}")


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
        print(self.aligned_sequences)
        for i, record in enumerate(self.aligned_sequences):
            if record.id == self.reference_sequence.id:
                self.reference_index = i
                break

    def set_CDS(self):
        CDS_dict = {}
        for feature in self.reference_sequence.features:
            if feature.type == 'CDS':
                for part in feature.location.parts:
                    gene = feature.qualifiers["gene"][0]
                    if gene not in CDS_dict:
                        CDS_dict[gene] = []
                    
                    # Check if the start position is already in the list(수정 가능성 있음)
                    for start, end in CDS_dict[gene]:
                        if int(part.start) == start:
                            break
                    # If not, add the start and end positions to the list
                    else:
                        CDS_dict[gene].append((int(part.start), int(part.end)))
        
        self.CDS_dict = CDS_dict
        
        for gene in CDS_dict:
            print(gene)
            for start, end in CDS_dict[gene]:
                print(start, end)


    def translate_sequences(self):
        for record in self.aligned_sequences:
            clean_dna_seq = str(record.seq).replace('-', 'N')
            clean_dna_seq = Seq(clean_dna_seq)
            protein_seq = clean_dna_seq.translate()
            self.protein_sequences[record.id] = SeqRecord(protein_seq, id=record.id)

    def write_protein_sequences(self):
        with open(self.protein_file, "w") as f:
            SeqIO.write(self.protein_sequences.values(), f, "fasta")

    def write_mutation(self):
        for record in self.variant_sequences:
            if self.reference_index is not None:
                mutation = [
                    i for i, (a, b) in enumerate(zip(self.protein_sequences[self.reference_sequence.id].seq, self.protein_sequences[record.id].seq))
                    if a != b and b != 'X'
                ]
                self.mutations.append((record.id, mutation))
            else:
                print("Reference sequence not found in alignment.")

    def get_mutation(self, file_handle):
        file_handle.write("Mutations:\n")
        for variant_id, mutation in self.mutations:
            file_handle.write(f"Mutations for variant {variant_id}:\n")
            file_handle.write(f"Total mutations for variant {variant_id}: {len(mutation)}/{len(self.protein_sequences[self.reference_sequence.id].seq)}\n")
            for mut in mutation:
                file_handle.write(f"Position {mut + 1}: {self.protein_sequences[self.reference_sequence.id].seq[mut]} -> {self.protein_sequences[variant_id].seq[mut]}\n")

    def get_metadata(self, file_handle):
        file_handle.write("Metadata:\n")
        file_handle.write(f"Sequence ID: {self.reference_sequence.id}\n")
        file_handle.write(f"Name: {self.reference_sequence.name}\n")
        file_handle.write(f"Description: {self.reference_sequence.description}\n")
        file_handle.write(f"Length: {len(self.reference_sequence)}\n")

    def get_gene_annotation(self, file_handle):
        file_handle.write("Gene annotations:\n")
        for feature in self.reference_sequence.features:
            if feature.type == "5'UTR" or feature.type == "3'UTR":
                file_handle.write(f"{feature.type}\n")
                file_handle.write(f"Location: {feature.location}\n\n")
            if feature.type == 'gene':
                file_handle.write(f"Gene: {feature.qualifiers['gene']}\n")
                file_handle.write(f"Location: {feature.location}\n\n")

    def get_peptide_annotation(self, file_handle):
        file_handle.write("Peptide annotations:\n")
        for feature in self.reference_sequence.features:
            if feature.type == 'mat_peptide':
                file_handle.write(f"Protein: {feature.qualifiers['product']}\n")
                file_handle.write(f"Location: {feature.location}\n\n")

    def get_CDS_annotation(self, file_handle):
        file_handle.write("CDS annotations:\n")
        for feature in self.reference_sequence.features:
            # print(feature)
            if feature.type == 'CDS':
                file_handle.write(f"Location: {feature}\n")
            # if feature.type == 'mat_peptide':
            #     file_handle.write(f"Protein: {feature.qualifiers['product']}\n")
            #     file_handle.write(f"Location: {feature.location}\n\n")

    def run(self):
        self.read_sequences()
        self.set_CDS()
        # self.run_muscle_dna()
        # self.read_alignment()
        # self.translate_sequences()
        # self.write_protein_sequences()
        # self.write_mutation()

        # with open("output.txt", "w") as file_handle:
            # self.get_metadata(file_handle)
            # self.get_gene_annotation(file_handle)
            # self.get_peptide_annotation(file_handle)
            # self.get_mutation(file_handle)
            # self.get_CDS_annotation(file_handle)

if __name__ == "__main__":
    reference_id = "NC_045512"
    files = ["MT576556.1.fasta"]
    alignment = SequenceAlignment(files, reference_id)
    alignment.run()
