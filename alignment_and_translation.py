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
        self.combined_file = "result/combined.fasta"
        self.aligned_file = "result/aligned.fasta"
        self.protein_file = "result/proteins.fasta"
        self.reference_sequence = None
        self.variant_sequences = []
        self.aligned_sequences = []
        self.protein_sequences = []
        self.mutations = []
        self.reference_index = None
        self.reference_id = reference_id
        self.CDS_dict = {}
        self.protein_dict={}
        self.alignment_dict={}

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

        except HTTPError as e:
            print(f"HTTPError: {e.code} - {e.reason}")


    def run_muscle_dna(self):
        result = subprocess.run([self.muscle_exe, "-in", self.combined_file, "-out", self.aligned_file])
        if result.returncode != 0:
            print("Error running muscle:")
            # print(result.stderr)
        else:
            print("Muscle ran successfully:")
            # print(result.stdout)

    def read_alignment(self):
        alignment = AlignIO.read(self.aligned_file, "fasta")
        desired_order = [self.reference_sequence.id] + [record.id for record in self.variant_sequences]
        self.alignment_dict = {record.id: record for record in alignment}
        self.aligned_sequences = [self.alignment_dict[id] for id in desired_order]


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
        
        # for gene in CDS_dict:
        #     print(gene)
        #     for start, end in CDS_dict[gene]:
        #         print(start, end)

    def translate_sequences(self):
        for gene, locations in  self.CDS_dict.items():
            for record in self.aligned_sequences:
                protein_seq = Seq("")
                for location in locations:
                    start = location[0]
                    end = location[1]
                    region_seq = record.seq[start:end]

                    # 예외처리 해야함 ex) Bio.Data.CodonTable.TranslationError: Codon '--G' is invalid 
                    protein_seq = protein_seq+region_seq.translate(to_stop=True)

                    # print(type(protein_seq))
                    # print(protein_seq)

                protein_record = SeqRecord(protein_seq, id=record.id, description=f"translated protein {start}-{end}")
                
                if gene not in self.protein_dict:
                    self.protein_dict[gene] = []
                self.protein_dict[gene].append(protein_record)

            with open(f"result/proteins_{gene}.fasta", "w") as f:
                SeqIO.write(self.protein_dict[gene], f, "fasta")

    # 재작성 필요
    # def set_mutation(self):
    #     for record in self.variant_sequences:
    #         if self.reference_index is not None:
    #             mutation = [
    #                 i for i, (a, b) in enumerate(zip(self.protein_sequences[self.reference_sequence.id].seq, self.protein_sequences[record.id].seq))
    #                 if a != b and b != 'X'
    #             ]
    #             self.mutations.append((record.id, mutation))
    #         else:
    #             print("Reference sequence not found in alignment.")

    def get_metadata(self):
        metadata = {
            "Sequence ID": self.reference_sequence.id,
            "Name": self.reference_sequence.name,
            "Description": self.reference_sequence.description,
            "Length": len(self.reference_sequence)
        }
        return metadata

    def get_gene_annotation(self):
        gene_annotations = []
        peptide_annotations = []

        for feature in self.reference_sequence.features:
            if feature.type == "5'UTR" or feature.type == "3'UTR":
                gene_annotations.append({
                    "type": feature.type,
                    "location": str(feature.location)
                })
            if feature.type == 'gene':
                gene_annotations.append({
                    "gene": feature.qualifiers['gene'],
                    "location": str(feature.location)
                })
            if feature.type == 'mat_peptide':
                peptide_annotations.append({
                    "protein": feature.qualifiers['product'],
                    "location": str(feature.location)
                })

        return {
            "Gene Annotations": gene_annotations,
            "Peptide Annotations": peptide_annotations
        }

    def run(self):
        self.read_sequences()
        self.set_CDS()
        self.run_muscle_dna()
        self.read_alignment()
        self.translate_sequences()
        # self.set_mutation()

if __name__ == "__main__":
    reference_id = "NC_045512"
    # files = ["data/MT576556.1.spike.fasta"]
    # files = ["data/MT576556.1.spike.fasta", "data/OR240434.1.spike.fasta", "data/PP346415.1.spike.fasta"]
    files = ["data/MT576556.1.fasta", "data/OR240434.1.fasta", "data/PP346415.1.fasta"]
    alignment = SequenceAlignment(files, reference_id)
    alignment.run()

    metadata = alignment.get_metadata()

    gene_annotation = alignment.get_gene_annotation()
    print(metadata)
    print(gene_annotation['Gene Annotations'])
