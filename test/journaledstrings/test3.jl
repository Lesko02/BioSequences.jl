include("journaledstrings.jl")
using BenchmarkTools

js1 = JournaledString(LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG"),
[SortedDict{Int, JournalEntry}() for _ in 1:10], 0)

add_delta!(js1, [1, 2], DeltaTypeIns, 8, "CGTA")
add_delta!(js1, [10], DeltaTypeSnp, 21, 'C')

println(js1.reference)
print_sequences(js1)
pattern = LongDNA{4}("CCC")
results1= Dict{Int64, Vector{UnitRange{Int64}}}()
results2 = Dict{Int64, Vector{UnitRange{Int64}}}()

#@benchmark slow_search(js1, pattern)
@benchmark approximate_search(js1, pattern)
#= println("NAIVE FIND")
print_results(results1)
println("EURISTICH FIND")
print_results(results2)=#