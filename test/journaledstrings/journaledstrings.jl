using BioSequences
using DataStructures
using FASTX

# Definition of DeltaTypes
@enum DeltaType DeltaTypeDel DeltaTypeIns DeltaTypeSnp DeltaTypeSV DeltaTypeCNV

# Custom insert for sequences
function insert!(seq::LongDNA{4}, pos::Int, subseq::LongDNA)
    for symbol in subseq
        Base.insert!(seq, pos, symbol)  # Inserimento dei singoli simboli
    end
    return seq
end

# Custom delete_at for ranges
function delete_at!(seq::LongDNA, pos_range::UnitRange{Int})
    for i in reverse(pos_range)
        Base.deleteat!(seq, i)
    end
    return seq
end

# Naive implementation of structure_variation //TO FIX//
function structure_variation!(seq::LongDNA{4}, pos::Int, subseq::LongDNA)
    # Check if it's a base case of an append!
    if pos == lastindex(seq)
        append!(seq, subseq)
        return seq
    # Check if it's a base case of insert!
    elseif pos < lastindex(seq)
        insert!(seq, pos, subseq)
        return seq
    # Check of an out of bounds insert --> SV
    elseif pos > lastindex(seq)
        i = lastindex(seq)
        range = lastindex(seq):pos-1
        for i in range
            push!(seq, '-')
        end
        append!(seq, subseq)
        return seq
    end
end

# Implementation of CNV
function copy_number_variation!(seq::LongDNA{4}, pos::Int,
     params::Tuple{LongDNA{4}, Int})
    subseq, rep = params  # Deconstruct the tuple
    for i in 1:rep  # Insert subseq `rep` times
        insert!(seq, pos, subseq)
    end    
    return seq
end


# Journal Entry
struct JournalEntry
    delta_type::DeltaType    # Type of delta
    position::Int            # Index of delta
    data::Any                # Additional Data
    time::Int                # Timestamp
end                          

# Definition of JournaledString Structure
mutable struct JournaledString
    reference::LongDNA{4}
    deltaMap::Vector{SortedDict{Int, JournalEntry, Base.Order.ForwardOrdering}}
    current_time::Int
end

const DeltaMap = SortedDict{Int, JournalEntry, Base.Order.ForwardOrdering}

# Constructor for JournaledString
function JournaledString(reference::LongDNA{4},
    deltaMap::Vector{DeltaMap})
    JournaledString(reference, deltaMap, 0)
end

struct JSTNode
    parent::Union{Nothing, JSTNode}
    deltaMap::Union{Nothing, DeltaMap}
    name::String
end

struct JSTree
    root::LongDNA{4};
    children::Dict{String, JSTNode}
 end

 function JSTree(root_sequence::LongDNA{4})
    # Create the root node with "Nothing" as the parent and an empty deltaMap
    root_node = JSTNode(nothing, nothing, "root")
    return JSTree(root_sequence, Dict("root" => root_node))
end



function add_node(tree::JSTree, parent_name::String, 
    deltas::DeltaMap, node_name::String)

    if !haskey(tree.children, parent_name)
        error("Parent node '$parent_name' does not exist.")
    else
        parent_node = tree.children[parent_name]
        new_node = JSTNode(parent_node, deltas, node_name)
    end

tree.children[node_name] = new_node
end


function flatten(tree::JSTree, node_name::String)
    # Base case
    if node_name == "root"
        return tree.root
    end

    # Recursive case
    node = tree.children[node_name]
    parent_sequence = flatten(tree, node.parent.name)
    return apply_delta(parent_sequence, node.deltaMap)
end

function print_tree(tree::JSTree, node_name::String = "root", indent::Int = 0)
    println(" "^(indent * 2) * "|- " * node_name) 
    for (child_name, child_node) in tree.children
        if child_node.parent !== nothing && child_node.parent.name == node_name
            print_tree(tree, child_name, indent + 1) 
        end
    end
end



function print_sequences(tree::JSTree)
    println("root: ", tree.root)
    for (name, node) in tree.children
        if name != "root"
        println("$name: ")
        println(flatten(tree, node.name))
        end
    end
end

function add_delta!(js::JournaledString, indices::Vector{Int}, 
                    delta_type::DeltaType, position::Int, data::Any)
    for idx in indices

       new_entry = JournalEntry(delta_type, position, data, js.current_time)
        
       js.deltaMap[idx][js.current_time] = new_entry

       js.current_time += 1
    end
end

function add_delta!(js::JournaledString,
     indices::Vector{Int}, entry::JournalEntry)

        for idx in indices

        js.deltaMap[idx][js.current_time] = entry

        js.current_time += 1
        end
end

function remove_delta!(js::JournaledString, time::Int)
    for idx in deltaMap
        for entry in js.deltaMap[idx]
            if entry[time] == time
                delete!(js.deltaMap[idx], time)
            else
            error("No mutation found at time: $time" )
            end
        end
    end
    
end

function apply_delta(reference::LongDNA{4}, delta::DeltaMap)
    seq = copy(reference)
    for (_, entry) in delta
        # Check DeltaType
        if entry.delta_type == DeltaTypeDel
            seq = delete_at!(seq, entry.position:(entry.position + 
                   entry.data - 1))  # Data is the bound of the range
        elseif entry.delta_type == DeltaTypeIns
            seq = insert!(seq, entry.position, LongDNA{4}(entry.data))
             # Single nucleotide permutation
        elseif entry.delta_type == DeltaTypeSnp
            seq[entry.position] = convert(DNA, entry.data)  
            # Larger Structure change
        elseif entry.delta_type == DeltaTypeSV
            seq = structure_variation!(seq, entry.position, entry.data)
        elseif entry.delta_type == DeltaTypeCNV
            seq = copy_number_variation!(seq, entry.position, entry.data)
        end
    end
    return seq
end

function print_sequences(jss::JournaledString)
    for i in 1:length(jss.deltaMap)
        modified_seq = apply_delta(jss.reference, jss.deltaMap[i])
        println("Sequence $i: ", modified_seq)
    end
end

# Function to build the n sequences into a sigle string
function build_sequences(jss::JournaledString)
    builded = ""
    for i in 1:length(jss.deltaMap)
        modified_seq = apply_delta(jss.reference, jss.deltaMap[i])
        builded *= "Sequence $i: " * string(modified_seq) * "\n"
    end
    return builded
end

# Function to print all of the deltas
function print_deltas(jss::JournaledString)
    for j in 1:length(jss.deltaMap)
        println("Stringa indice $j:")
        for (time, entry) in jss.deltaMap[j]
            println("  [time=$time] JournalEntry: $entry")
        end
    end
end

#API definitions

function get_mutation_history(jss::JournaledString)
    mutation_history = ""
    for (time, entry) in jss.delta_map
        mutation_history *= "Time $time: $entry\n"
    end
    return mutation_history
end

function get_mutation_interval(jss::JournaledString, time1::Int, time2::Int)
    mutation_interval = ""
    for (time, entry) in jss.delta_map
        if time1 <= time <= time2
            mutation_interval *= "Time $time: $entry\n"
        end
    end
    return mutation_interval
end

function get_sequences_at_time(jss::JournaledString, time::Int)
    sequences_at_time = Vector{LongDNA{4}}(undef, length(jss.deltaMap))
    for i in 1:length(jss.deltaMap)
        filtered_delta = DeltaMap()
        for (entry_time, entry) in jss.deltaMap[i]
            if entry_time > time
                break
            end
            filtered_delta[entry_time] = entry
        end
        sequences_at_time[i] = apply_delta(jss.reference, filtered_delta)
    end
    return sequences_at_time
end

function simulate_mutation!(jss::JournaledString, index::Int,
                             entry::JournalEntry)
   add_delta!(jss.deltaMap, index, entry)
end

function remove_mutation!(jss::JournaledString, time::Int)
    remove_delta!(jss.deltaMap, time)
end

function is_equal(jst1::JournaledString, jst2::JournaledString)::Bool
    if hash(jst1.reference) != hash(jst2.reference)
        return false
    end
    if jst1.current_time != jst2.current_time
        return false
    end
return hash(jst1.deltaMap)==hash(jst2.deltaMap)
end

function exact_search(jss::JournaledString, needle::LongDNA )
    results = Dict(i => UnitRange{Int64}[] for i in 1:length(jss.deltaMap))
    query = ExactSearchQuery(needle)
    vector = UnitRange{Int64}[]
    for i in 1:length(jss.deltaMap)

        seq = apply_delta(jss.reference, jss.deltaMap[i])
        vector = BioSequences.findall(query, seq)
        append!(results[i], vector)

    end
    return results
end

function apply_delta(reference::LongDNA{4}, entry::JournalEntry)
    seq = copy(reference)
        # Check on the DeltaType
        if entry.delta_type == DeltaTypeDel
            seq = delete_at!(seq, entry.position:(entry.position + 
            entry.data - 1))  # Data is the bound of the range
        elseif entry.delta_type == DeltaTypeIns
            seq = insert!(seq, entry.position, LongDNA{4}(entry.data))
        # Single nucleotide permutation
        elseif entry.delta_type == DeltaTypeSnp
            seq[entry.position] = convert(DNA, entry.data)  
            # Larger Structure change
        elseif entry.delta_type == DeltaTypeSV
            seq = structure_variation!(seq, entry.position, entry.data)
        elseif entry.delta_type == DeltaTypeCNV
            seq = copy_number_variation!(seq, entry.position, entry.data)
        end
    return seq
end

function exact_search(jst::JSTree, needle::LongDNA{4})
    indexMatrix = Dict{String, Vector{UnitRange{Int64}}}()
    query = ExactSearchQuery(needle)
    vector = UnitRange{Int}[]

    indexMatrix = Dict{String, Vector{UnitRange{Int64}}}(
        name => UnitRange{Int64}[] for name in keys(jst.children))
    
    for (name, child) in jst.children

        empty!(vector)
        if (!isnothing(child.deltaMap))
        seq = flatten(jst, name)
        vector = BioSequences.findall(query, seq)
        end

        indexMatrix[name] = append!(indexMatrix[name], vector)
        
    end
    return indexMatrix
end

function approximate_findall(query, tolerance::Int64, seq::LongDNA{4})
    results = UnitRange{Int64}[]
    pos = findfirst(query, tolerance, seq)
    while pos !== nothing
        push!(results, pos)
        pos = findnext(query, tolerance, seq, last(pos)+tolerance+1)
    end
    return results
end

function approximate_search(jss::JournaledString, needle::LongDNA{4})
    query = ApproximateSearchQuery(needle)
    tolerance = ceil(Int64, (length(needle) / 100) * 5)  
    indexMatrix = Dict{Int64, Vector{UnitRange{Int64}}}()
    vector = approximate_findall(query, tolerance, jss.reference)
    to_remove = Set{UnitRange{Int64}}()
    to_add = Set{UnitRange{Int64}}()

    for i in 1:length(jss.deltaMap)
        indexMatrix[i] = vector
    end

    for i in 1:length(jss.deltaMap)
        empty!(to_add)
        empty!(to_remove)
        for range in indexMatrix[i]
            for ( _, entry) in jss.deltaMap[i]
                
                if entry.position in range
                    seq = apply_delta(jss.reference, entry)

                    for element in approximate_findall(query, tolerance, seq)
                        push!(to_add, element)
                    end
                    
                    push!(to_remove, range)
                end

            end
        end 
        indexMatrix[i]= filter(x -> all(y -> x != y, to_remove), indexMatrix[i])
        append!(indexMatrix[i], to_add)
        indexMatrix[i] = collect(Set(indexMatrix[i]))
    end
    return indexMatrix
end

function print_results(results::Dict{Int64, Vector{UnitRange{Int64}}})
    for i in 1:length(results)
        if isempty(results[i])
            println("No Match in DeltaMap $i")
        else
            println("Match in $i:")
            println("Ranges: ", results[i])
        end
    end
end

function print_results(results::Dict{String, Vector{UnitRange{Int64}}})
    for (name, _) in results
        if isempty(results[name])
            println("No Match at $name")
        else
            println("Match at $name:")
            println("Ranges: ", results[name])
        end
    end
end

function approximate_search(jst::JSTree, needle::LongDNA{4})
    query = ApproximateSearchQuery(needle)
    indexes = UnitRange{Int}[]
    vector = UnitRange{Int}[]
    to_remove = UnitRange{Int64}[]   
    indexes = findall(query, seq)
    tolerance = ceil(Int64, (length(needle) / 100) * 5) 

   for range in indexes
        for (name, node) in jst.children
            for ( _, entry) in node.deltaMap
                if entry.position in range
                    seq = flatten(jst, name)
                    empty!(vector)
                    seq = apply_delta(seq, entry)
                    vector = push!(approximate_findall(query, tolerance, seq))
                    push!(to_remove, range)
                    if isempty(vector)
                        println("No match at Child $name")
                    else
                        println("Match at Child $name")
                        println("Ranges: ", vector)
                    end
                end
            end
            filter!(x -> x ∉ to_remove, indices)
            if !isempty(indices)
                println("Match at Child $name")
                println("Ranges: ", indices) 
            end
        end
    end   
end

function approximate_search(jss::JournaledString, needle::LongDNA{4},
    tol::Int64)

if tol <= 0 || tol >= 100
    error("Tolerance cannot less or 0% or more than 100%")
end

tolerance = ceil(Int64, (length(needle) / 100) * tol) 
query = ApproximateSearchQuery(needle)
indexMatrix = Dict{Int64, Vector{UnitRange{Int64}}}()
vector = approximate_findall(query, tolerance, jss.reference)
to_remove = Set{UnitRange{Int64}}()
to_add = Set{UnitRange{Int64}}()

for i in 1:length(jss.deltaMap)
    indexMatrix[i] = vector
end

for i in 1:length(jss.deltaMap)
    empty!(to_add)
    empty!(to_remove)
    for range in indexMatrix[i]
        for ( _, entry) in jss.deltaMap[i]
            
            if entry.position in range
                seq = apply_delta(jss.reference, entry)

                for element in approximate_findall(query, tolerance, seq)
                    push!(to_add, element)
                end
                
                push!(to_remove, range)
            end

        end
    end 
    indexMatrix[i]= filter(x -> all(y -> x != y, to_remove), indexMatrix[i])
    append!(indexMatrix[i], to_add)
    indexMatrix[i] = collect(Set(indexMatrix[i]))
end
return indexMatrix
end