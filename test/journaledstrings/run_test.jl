include("journaledstrings.jl")
# Define a reference sequence
reference_seq = LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG")

# Initialize the deltaMap for how many sequences you want
deltaMap = [SortedDict{Int, JournalEntry}() for _ in 1:10]  # Ten subseqs

# Create a JournaledString
jst = JournaledString(reference_seq, deltaMap)

println(typeof(jst))
# Example: Apply insertions, deletions, SNPs, etc.
println("Prova di: INS e SNP")
add_delta!(jst, [1, 2], DeltaTypeIns, 8, "CGTA")
add_delta!(jst, [4, 9], DeltaTypeSnp, 10, 'C')
add_delta!(jst, [8], DeltaTypeIns, 24, "NNNNN")
# print_sequences(jst)

println("\n")
println("Prova di: DELETE (stringa 3 perde 10 caratteri)")
add_delta!(jst, [3], DeltaTypeDel, 1, 10)
# print_sequences(jst)

println("\n")
println("Prova di: structure_variation")
add_delta!(jst, [5], DeltaTypeSV, 30, dna"CGTACGTACGTACGTA")
add_delta!(jst, [10], DeltaTypeSV, 18, dna"CGTACGTACGTACGTA")
add_delta!(jst, [7], DeltaTypeSV, 24, dna"NNNNN")
add_delta!(jst, [7], DeltaTypeSV, 24, dna"NNNNN")
add_delta!(jst, [7], DeltaTypeSV, 24, dna"NNNNN")
add_delta!(jst, [7], DeltaTypeSV, 24, dna"NNNNN")
add_delta!(jst, [7], DeltaTypeSV, 24, dna"NNNNN")

# Print the sequences and the deltas
# print_sequences(jst)

println("\n")
println("Prova di: CNV")
add_delta!(jst, [1], DeltaTypeCNV, 1, (LongDNA{4}("AAAA"), 5))

print_sequences(jst)
print_deltas(jst)

# Create two identical JournaledString objects
js1 = JournaledString(LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG"), [SortedDict{Int, JournalEntry}() for _ in 1:10], 0)
js2 = JournaledString(LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG"), [SortedDict{Int, JournalEntry}() for _ in 1:10], 0)
add_delta!(js1, [1, 2], DeltaTypeIns, 8, "CGTA")
add_delta!(js2, [1, 2], DeltaTypeIns, 8, "CGTA")
# Compare two JournaledString objects
println(is_equal(js1, js2))


tree = JSTree(LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG"))

add_node(tree, "root", deltaMap[2], "child1")
add_node(tree, "child1", deltaMap[1], "child2")
add_node(tree, "child1", deltaMap[3], "child3")
add_node(tree, "child1", deltaMap[4], "child4")
add_node(tree, "child4", deltaMap[5], "child5")

print_tree(tree, "root", 0)
print_sequences(tree)