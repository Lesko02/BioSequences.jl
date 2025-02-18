include("journaledstrings.jl")
using BenchmarkTools
# Define a reference sequence
reference_seq = LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG")

deltaMap = [SortedDict{Int, JournalEntry}() for _ in 1:5]  # Ten subseqs

# Create a JournaledString
jst = JournaledString(reference_seq, deltaMap)

add_delta!(jst, [1, 2], DeltaTypeIns, 8, "CGTA")
add_delta!(jst, [4], DeltaTypeSnp, 10, 'C')
add_delta!(jst, [3], DeltaTypeDel, 1, 10)
add_delta!(jst, [5], DeltaTypeSV, 30, dna"CGTACGTACGTACGTACCC")
add_delta!(jst, [4], DeltaTypeSV, 24, dna"NNNNN")
add_delta!(jst, [1], DeltaTypeCNV, 1, (LongDNA{4}("ATCG"), 2))



# Create two identical JournaledString objects
js1 = JournaledString(LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG"),
 deltaMap, 0)

tree = JSTree(LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG"))
add_node(tree, "root", deltaMap[1], "child1")
add_node(tree, "child1", deltaMap[2] , "child2")
add_node(tree, "child1", deltaMap[3], "child3")
add_node(tree, "child1", deltaMap[4], "child4")
add_node(tree, "child4", deltaMap[5], "child5")

################################
#print_sequences(js1)
pattern = LongDNA{4}("CCC")
# @time exact_search(js1, pattern)
################################
println("\n")
print_tree(tree)
print_sequences(tree)
results = approximate_search(tree, pattern)
println("\n")
print_results(results)


