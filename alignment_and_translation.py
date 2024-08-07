from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, Entrez
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from urllib.error import HTTPError
import subprocess
import os


class SequenceAlignment:
    def __init__(self, files, reference_id, muscle_exe="muscle"):
        Entrez.email = "your_email@example.com"
        self.files = files
        self.muscle_exe = muscle_exe
        self.combined_file = "result/combined.fasta"
        self.aligned_file = "result/aligned.fasta"
        self.reference_sequence = None
        self.variant_sequences = []
        self.aligned_sequences = []
        self.reference_id = reference_id
        self.alignment_dict={}
        self.protein_length={}
        self.mutation_dict={}
        self.reference_protein_seq=""
        self.alignment_index={}
        
        try:
            # Get reference sequence
            handle = Entrez.efetch(db="nucleotide", id=reference_id, rettype="gb", retmode="text")
            self.reference_sequence = SeqIO.read(handle, "genbank")
            handle.close()
            self.metadata={
                "Sequence ID": self.reference_sequence.id,
                "Name": self.reference_sequence.name,
                "Description": self.reference_sequence.description,
                "Length": len(self.reference_sequence)
            }
        except HTTPError as e:
            print(f"HTTPError: {e.code} - {e.reason}")
            

    def read_sequences(self):
        # Get reference protein sequence
        CDS_dict={}
        for feature in self.reference_sequence.features:
            if feature.type == 'CDS':
                gene = feature.qualifiers["gene"][0]
                if gene not in CDS_dict:
                    CDS_dict[gene] = feature.qualifiers["translation"][0]
                
        for gene, protein_seq in CDS_dict.items():
            self.protein_length[gene] = len(protein_seq)
            self.reference_protein_seq += protein_seq

        # Get variant sequences
        self.variant_sequences = [SeqIO.read(file, "fasta") for file in self.files]
        translated_variants = [variant.seq.translate(to_stop=True) for variant in self.variant_sequences]
        translated_records = [SeqRecord(Seq(translation), id=variant.id, description="translated variant protein")
                              for variant, translation in zip(self.variant_sequences, translated_variants)]
        
        # Combine reference and translated variant sequences
        reference_protein_record = SeqRecord(Seq(self.reference_protein_seq), id=self.reference_sequence.id, description="reference protein")
        with open(self.combined_file, "w") as f:
            SeqIO.write([reference_protein_record] + translated_records, f, "fasta")

    def run_muscle_dna(self):
        # Run muscle
        result = subprocess.run([self.muscle_exe, "-in", self.combined_file, "-out", self.aligned_file])
        if result.returncode != 0:
            print("Error running muscle:")
        else:
            print("Muscle ran successfully:")


    def read_alignment(self):
        # Read alignment
        alignment = list(SeqIO.parse(self.aligned_file, "fasta"))

        # Sort alignment
        desired_order = [self.reference_sequence.id] + [record.id for record in self.variant_sequences]
        self.alignment_dict = {record.id: record for record in alignment}
        self.aligned_sequences = [self.alignment_dict[id] for id in desired_order]

        # Update protein length
        record = self.aligned_sequences[0]
        start=0
        for gene, length in self.protein_length.items():
            end = start + length
            gap_count = record.seq[start:end].count('-')
            end += gap_count
            self.protein_length[gene] = length+gap_count

            # Set alignment data
            self.alignment_index[gene] = (start, end)

            start=end

        ########################################
        ###### Peptide 구간 코드 구현 필요 #######
        ### if feature.type == 'mat_peptide': ##
        ########################################

    def write_protein_sequences(self):
        # for (gene, start, end) in self.alignment_index:
        for gene, (start, end) in self.alignment_index.items():
            protein_sequences=[]
            for record in self.aligned_sequences:
                protein_record = SeqRecord(record.seq[start:end], id=record.id, description=f"translated protein {start}-{end}")
                protein_sequences.append(protein_record)

            with open(f"result/proteins_{gene}.fasta", "w") as f:
                SeqIO.write(protein_sequences, f, "fasta")


    def set_mutation(self):
        reference_protein = self.aligned_sequences[0].seq
        for record in self.aligned_sequences[1:]:
            variant_protein = record.seq
            mutation = []
            for i, (ref, var) in enumerate(zip(reference_protein, variant_protein)):
                if ref != var and ref != "-" and var != "-":
                    mutation.append((i, ref, var))
            self.mutation_dict[record.id] = mutation

    def run_linear_design(self, gene, variant_id):
        (start,end) = self.alignment_index[gene]
        input_sequence = str(self.alignment_dict[variant_id].seq[start:end]).replace("-", "")
        
        # for test 긴 길이 통으로 분석 안됨
        input_sequence=input_sequence[318:541]
        os.chdir("./LinearDesign")
        command = f"echo {input_sequence} | ./lineardesign"
        exit_code = os.system(command)

        if exit_code == 0:
            print("Command executed successfully")
        else:
            print("Error executing command")
        os.chdir("..")




    def get_metadata(self):
        return self.metadata
    
    def get_alignment_data(self):
        return self.alignment_index, self.aligned_sequences

    def get_mutation(self):
        return self.mutation_dict
    
    def get_prot_param(self):
        # 분석할 단백질 서열
        sequence = "RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF"

        # 단백질 서열 분석 객체 생성
        protein_analysis = ProteinAnalysis(sequence)

        # 분자량 계산
        molecular_weight = protein_analysis.molecular_weight()
        print(f"Molecular Weight: {molecular_weight:.2f} Da")

        # 아미노산의 개수
        amino_acid_count = protein_analysis.count_amino_acids()
        print("Amino Acid Count:")
        for aa, count in amino_acid_count.items():
            print(f"{aa}: {count}")

        # 아미노산 비율
        amino_acid_percent = protein_analysis.get_amino_acids_percent()
        print("Amino Acid Percent:")
        for aa, percent in amino_acid_percent.items():
            print(f"{aa}: {percent:.2%}")

        # 이성질화점 (pI)
        isoelectric_point = protein_analysis.isoelectric_point()
        print(f"Isoelectric Point (pI): {isoelectric_point:.2f}")

        # 불안정성 지수
        instability_index = protein_analysis.instability_index()
        print(f"Instability Index: {instability_index:.2f}")

        # 극성, 비극성, 기본성, 산성 아미노산 비율
        secondary_structure_fraction = protein_analysis.secondary_structure_fraction()
        print(f"Secondary Structure Fraction (Helix, Turn, Sheet): {secondary_structure_fraction}")

        # 향수성 지수
        gravy = protein_analysis.gravy()
        print(f"Gravy: {gravy:.2f}")

        # 극성 잔기 비율
        aromaticity = protein_analysis.aromaticity()
        print(f"Aromaticity: {aromaticity:.2%}")

    def run(self):
        self.read_sequences()
        self.run_muscle_dna()
        self.read_alignment()
        self.write_protein_sequences()    # 부가기능(LinearDesign) 활용 시 주석 해제
        self.set_mutation()
        # self.run_linear_design("S", "MW642250.1")

if __name__ == "__main__":
    reference_id = "NC_045512"
    files = ["data/OL672836.1.spike.fasta", "data/MW642250.1.spike.fasta", "data/OM958567.1.spike.fasta"]
    alignment = SequenceAlignment(files, reference_id)

    # alignment 실행 전 metadata 받아오기
    metadata = alignment.get_metadata()

    # alignment 실행
    alignment.run()

    # alignment data 받아오기
    alignment_index, aligned_sequences = alignment.get_alignment_data()
  
    # mutation data 받아오기
    mutation_dict = alignment.get_mutation()

    ########################################
    ###### metadata, alignment data 확인 ####
    ########################################
    # Example for use data
    # for gene, (start, end) in alignment_index.items():
    #     print(gene, start, end)
    # for record in aligned_sequences:
    #     print(record.id)
    #     print(record.seq)
    #     print()
    # for key, value in mutation_dict.items():
    #     print(key)
    #     for i, ref, var in value:
    #         print(f"{i}: {ref} -> {var}")
    #     print()
