from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO, SeqIO, Entrez
from urllib.error import HTTPError
import subprocess
from Bio import Align


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
        # self.protein_sequences = []
        self.mutations = []
        self.reference_index = None
        self.reference_id = reference_id
        self.CDS_dict = {}
        self.alignment_dict={}
        self.protein_dict={}
        self.mutation_dict={}
        self.concatenated_protein_seq=""
        self.len_list = []
        
        try:
            # Get reference sequence
            handle = Entrez.efetch(db="nucleotide", id=reference_id, rettype="gb", retmode="text")
            self.reference_sequence = SeqIO.read(handle, "genbank")
            handle.close()
        except HTTPError as e:
            print(f"HTTPError: {e.code} - {e.reason}")

    def set_reference_protein(self):
        for feature in self.reference_sequence.features:
            if feature.type == 'CDS':
                gene = feature.qualifiers["gene"][0]
                if gene not in self.CDS_dict:
                    self.CDS_dict[gene] = feature.qualifiers["translation"][0]
                
        for gene, protein_seq in self.CDS_dict.items():
            # protein_record = SeqRecord(Seq(protein_seq[0]), id=self.reference_sequence.id, description=f"reference protein")
            
            self.protein_dict[gene] = [protein_seq, len(protein_seq)]
            self.concatenated_protein_seq += protein_seq
        protein_record = SeqRecord(Seq(self.concatenated_protein_seq), id=self.reference_sequence.id, description="reference protein")
        
        with open("concatenated_protein.fasta", "w") as f:
            SeqIO.write(protein_record, f, "fasta")
        
        # for gene, protein_seq in self.protein_dict.items():
        #     print(f"Gene: {gene}, Protein: {protein_seq[0]}, Length: {protein_seq[1]}")


    def read_sequences(self):
        # Get variant sequences
        self.variant_sequences = [SeqIO.read(file, "fasta") for file in self.files]
        translated_variants = [variant.seq.translate(to_stop=True) for variant in self.variant_sequences]
        translated_records = [SeqRecord(Seq(translation), id=variant.id, description="translated variant protein")
                              for variant, translation in zip(self.variant_sequences, translated_variants)]
        
        # Combine reference and translated variant sequences
        # test Spike only
        # self.concatenated_protein_seq = self.CDS_dict["S"]
        reference_protein_record = SeqRecord(Seq(self.concatenated_protein_seq), id=self.reference_sequence.id, description="reference protein")
        with open(self.combined_file, "w") as f:
            SeqIO.write([reference_protein_record] + translated_records, f, "fasta")

        # 확인을 위해 출력
        # for record in translated_records:
        #     print(record.format("fasta"))



    def run_muscle_dna(self):
        result = subprocess.run([self.muscle_exe, "-in", self.combined_file, "-out", self.aligned_file])
        if result.returncode != 0:
            print("Error running muscle:")
        else:
            print("Muscle ran successfully:")


    def read_alignment(self):
        alignment = AlignIO.read(self.aligned_file, "fasta")
        desired_order = [self.reference_sequence.id] + [record.id for record in self.variant_sequences]
        self.alignment_dict = {record.id: record for record in alignment}
        self.aligned_sequences = [self.alignment_dict[id] for id in desired_order]



    def set_CDS(self):
        CDS_dict = {}
        for feature in self.reference_sequence.features:
            # Check if feature is CDS and handle appropriately
            if feature.type == 'CDS':
                gene = feature.qualifiers["gene"][0]
                if gene not in CDS_dict:
                    CDS_dict[gene] = {
                        'CDS': [],
                        'sequence': '',
                        # 'peptide': []
                    }
                start = int(feature.location.start)
                # Check for duplicate CDS start positions
                if any(start == cds[0] for cds in CDS_dict[gene]['CDS']):
                    continue  # Skip duplicate CDS
                CDS_dict[gene]['CDS'].append((start, int(feature.location.end)))
                CDS_dict[gene]['sequence'] += feature.qualifiers["translation"][0]
            # Check for peptides and add if valid (not after duplicate CDS)
            # elif feature.type == 'mat_peptide':
            #     # print(feature)
            #     gene = feature.qualifiers["gene"][0]
            #     # Only add peptide if it follows the most recent CDS without duplication
            #     peptide_info = (feature.qualifiers["product"][0], int(feature.location.start), int(feature.location.end))
            #     CDS_dict[gene]['peptide'].append(peptide_info)
        
        self.CDS_dict = CDS_dict
        
        # print(CDS_dict)
        # for gene in CDS_dict:
            # print(gene)
            # for start, end in CDS_dict[gene]:
            #     print(start, end)

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

    
    def set_mutation(self):
        for gene, records in self.protein_dict.items():
            reference_record = None
            mutation_record = {}
            for record in records:
                if(record.id == self.reference_sequence.id):
                    reference_record = record
                    continue
                mutations = []
                for i, (ref_aa, var_aa) in enumerate(zip(reference_record.seq, record.seq)):
                    if ref_aa != var_aa and var_aa != '-':
                        mutations.append((i, ref_aa, var_aa))
                
                mutation_record[record.id] = mutations
            self.mutation_dict[gene] = mutation_record
        
        # mutation_dict 출력
        for gene, mutations_record in self.mutation_dict.items():
            print(f"Gene: {gene}")
            for variant_id, mutations in mutations_record.items():
                if len(mutations) == 0:
                    continue
                print(f"Variant: {variant_id}")
                for mutation in mutations:
                    print(f"Position: {mutation[0]}, Reference: {mutation[1]}, Variant: {mutation[2]}")
                


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
                    "gene": feature.qualifiers['gene'][0],
                    "location": str(feature.location)
                })
            if feature.type == 'mat_peptide':
                peptide_annotations.append({
                    "protein": feature.qualifiers['product'][0],
                    "location": str(feature.location)
                })

        return {
            "Gene Annotations": gene_annotations,
            "Peptide Annotations": peptide_annotations
        }

    def run(self):
        self.set_reference_protein()
        self.read_sequences()
        self.set_CDS()
        # self.run_muscle_dna()
        # self.read_alignment()
        # self.translate_sequences()
        # self.set_mutation()

if __name__ == "__main__":
    reference_id = "NC_045512"
    # files = ["data/MT576556.1.spike.fasta"]
    # files = ["data/MT576556.1.spike.fasta", "data/OR240434.1.spike.fasta", "data/PP346415.1.spike.fasta"]
    # files = ["data/MT576556.1.fasta", "data/OR240434.1.fasta", "data/PP346415.1.fasta"]
    files = ["data/OL672836.1.spike.fasta", "data/MW642250.1.spike.fasta"]
    alignment = SequenceAlignment(files, reference_id)
    alignment.run()

    # metadata = alignment.get_metadata()

    # gene_annotation = alignment.get_gene_annotation()
    # print(metadata)
    # print(gene_annotation)
    # for gene in gene_annotation['Gene Annotations']:
        # print(gene)
