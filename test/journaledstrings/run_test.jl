include("journaledstrings.jl")
# Define a reference sequence
reference_seq = LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG")

# Initialize the deltaMap for how many sequences you want
deltaMap = [SortedDict{Int, JournalEntry}() for _ in 1:10]  # Ten subseqs


# Create a JournaledString
jst = JournaledString(reference_seq, deltaMap)

# Example: Apply insertions, deletions, SNPs, etc.
println("Prova di: INS e SNP")
add_delta!(jst.deltaMap, [1, 2], DeltaTypeIns, 8, "CGTA")
add_delta!(jst.deltaMap, [4, 9], DeltaTypeSnp, 10, 'C')
add_delta!(jst.deltaMap, [8], DeltaTypeIns, 24, "NNNNN")
print_sequences(jst)

println("\n")
println("Prova di: DELETE (stringa 3 perde 10 caratteri)")
add_delta!(jst.deltaMap, [3], DeltaTypeDel, 1, 10)
print_sequences(jst)

println("\n")
println("Prova di: structure_variation")
add_delta!(jst.deltaMap, [5], DeltaTypeSV, 30, dna"CGTACGTACGTACGTA")
add_delta!(jst.deltaMap, [10], DeltaTypeSV, 18, dna"CGTACGTACGTACGTA")
add_delta!(jst.deltaMap, [7], DeltaTypeSV, 24, dna"NNNNN")
add_delta!(jst.deltaMap, [7], DeltaTypeSV, 24, dna"NNNNN")
add_delta!(jst.deltaMap, [7], DeltaTypeSV, 24, dna"NNNNN")
add_delta!(jst.deltaMap, [7], DeltaTypeSV, 24, dna"NNNNN")
add_delta!(jst.deltaMap, [7], DeltaTypeSV, 24, dna"NNNNN")

# Print the sequences and the deltas
print_sequences(jst)
print_deltas(jst)
