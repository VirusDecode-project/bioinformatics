from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# 예제 DNA 서열 (reference와 variant)
reference_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
variant_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")

# DNA 서열을 단백질 서열로 번역
reference_protein = reference_dna.translate()
variant_protein = variant_dna.translate()

print("Reference Protein Sequence:")
print(reference_protein)
print("Variant Protein Sequence:")
print(variant_protein)

# 서열 정렬
alignments = pairwise2.align.globalxx(reference_protein, variant_protein)

# 가장 좋은 정렬 결과 출력
print("Best Alignment:")
print(format_alignment(*alignments[0]))

# 뮤테이션 찾기
ref_aligned = alignments[0][0]
var_aligned = alignments[0][1]

mutations = []
for i in range(len(ref_aligned)):
    if ref_aligned[i] != var_aligned[i]:
        mutations.append((i, ref_aligned[i], var_aligned[i]))

print("Mutations found:")
for (pos, ref, var) in mutations:
    print(f"Position: {pos}, Reference: {ref}, Variant: {var}")
