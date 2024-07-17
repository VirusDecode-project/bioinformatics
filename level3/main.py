from Bio import Entrez, SeqIO

# NCBI에 접근할 이메일 설정 (NCBI 정책에 따라 반드시 필요함)
Entrez.email = "your_email@example.com"

def fetch_annotation(sequence_id):
    try:
        # NCBI에서 염기서열 정보 가져오기
        handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        x=1
        # 염기서열 어노테이션 정보 출력
        print(f"Sequence ID: {record.id}")
        print(f"Name: {record.name}")
        print(f"Description: {record.description}")
        print("Annotations:")
        for feature in record.features:
            # print(feature)
            if feature.type == "5'UTR" or feature.type =="3'UTR":
                print(feature.type)
                print(f"Location: {feature.location}")
                print()
            if feature.type == 'gene':
                print(f"Gene: {feature.qualifiers['gene']}")
                print(f"Location: {feature.location}")
                print()
            # if feature.type == 'mat_peptide':
            #     print(f"Protein: {feature.qualifiers['product']}")
            #     print(f"Location: {feature.location}")
            #     print()

    except Exception as e:
        print(f"An error occurred: {e}")

# 특정 염기서열 ID (예: "NC_045512")를 이용하여 어노테이션 정보 가져오기
fetch_annotation("NC_045512.2")
