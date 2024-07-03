from Bio import Entrez, SeqIO

# NCBI에 이메일 주소 등록
Entrez.email = "example@gmail.com"

class NCBISequenceFetcher:
    def fetch_nucleotide_sequence(self, nucleotide_id):
        handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="gb", retmode="text")
        seq_record = SeqIO.read(handle, "genbank")
        handle.close()
        return seq_record

    def fetch_fasta_sequence(self, nucleotide_id, file_path):
        handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="fasta", retmode="text")
        with open(file_path, 'w') as fasta_file:
            fasta_file.write(handle.read())
        handle.close()
        print(f"FASTA sequence saved to {file_path}")

class SequenceMetadata:
    def __init__(self, seq_record):
        self.seq_record = seq_record

    def get_id(self):
        return self.seq_record.id

    def get_name(self):
        return self.seq_record.name

    def get_description(self):
        return self.seq_record.description

    def get_length(self):
        return len(self.seq_record)

    def get_sequence(self):
        return self.seq_record.seq

    def print_sequence_info(self):
        print(f"ID: {self.get_id()}")
        print(f"Name: {self.get_name()}")
        print(f"Description: {self.get_description()}")
        print(f"Length: {self.get_length()}")
    
    def print_sequence(self):
        print(f"Sequence: {self.get_sequence()}")

# 사용 예시
nucleotide_id = "NC_045512"  # 예: coronavirus 2

# 서열 데이터 가져오기
fetcher = NCBISequenceFetcher()
sequence_record = fetcher.fetch_nucleotide_sequence(nucleotide_id)

# 서열의 메타데이터 출력
metadata = SequenceMetadata(sequence_record)
metadata.print_sequence_info()

# FASTA 파일 저장
fetcher.fetch_fasta_sequence(nucleotide_id, "sequence.fasta")
