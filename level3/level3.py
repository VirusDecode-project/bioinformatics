from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# FASTA 파일 경로
fasta_file = "NC_045512.fasta"

# 염기서열 읽기
record = SeqIO.read(fasta_file, "fasta")
sequence = record.seq

print("1")
# NCBI BLAST 실행
result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
print("2")
# BLAST 결과 파싱
blast_records = NCBIXML.parse(result_handle)
print("3")
# BLAST 결과 분석
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print(f"****Alignment****")
            print(f"sequence: {alignment.title}")
            print(f"length: {alignment.length}")
            print(f"e value: {hsp.expect}")
            print(f"{hsp.query[0:75]}...")
            print(f"{hsp.match[0:75]}...")
            print(f"{hsp.sbjct[0:75]}...")
