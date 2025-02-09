include("journaledstrings.jl")
using BenchmarkTools

js1 = JournaledString(LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG"),
[SortedDict{Int, JournalEntry}() for _ in 1:10], 0)

add_delta!(js1, [1, 2], DeltaTypeIns, 8, "CGTA")
add_delta!(js1, [10], DeltaTypeSnp, 21, 'C')
add_delta!(js1, [4], DeltaTypeSV, 24, dna"NNNNN")
add_delta!(js1, [1, 3], DeltaTypeCNV, 1, (LongDNA{4}("ATCG"), 2))
add_delta!(js1, [5, 6, 7], DeltaTypeIns, 5, "GTC")
add_delta!(js1, [8, 9], DeltaTypeDel, 10, 2)
add_delta!(js1, [4, 2], DeltaTypeSnp, 18, 'T')
add_delta!(js1, [1, 3], DeltaTypeSnp, 18, 'G')
add_delta!(js1, [6, 8], DeltaTypeSV, 20, dna"CCTG")
add_delta!(js1, [5, 7, 10], DeltaTypeDel, 3, 2)

println(js1.reference)
print_sequences(js1)
pattern = LongDNA{4}("ATCG")
#prova= approximate_search(js1, pattern)
#print_results(prova)
#@benchmark slow_search(js1, pattern)
@benchmark approximate_search(js1,pattern)

