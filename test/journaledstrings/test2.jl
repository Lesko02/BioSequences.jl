include("journaledstrings.jl")
using BenchmarkTools

js1 = JournaledString(LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG"),
[SortedDict{Int, JournalEntry}() for _ in 1:10], 0)

add_delta!(js1, [1, 2], DeltaTypeIns, 8, "CGTA")
add_delta!(js1, [10], DeltaTypeSnp, 21, 'C')
pattern = LongDNA{4}("CCC")
@time slow_search(js1,pattern)
@time approximate_search(js1,pattern)
@btime slow_search(js1,pattern)
@btime approximate_search(js1,pattern)
