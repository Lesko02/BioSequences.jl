using BioSequences

# Custom insert per sequenze 
function insert!(seq::LongDNA{4}, pos::Int, subseq::LongDNA)
    for symbol in subseq
        Base.insert!(seq, pos, symbol)  # Inserimento dei singoli simboli
    end
    return seq
end

# Custom delete_at per intervalli 
function delete_at!(seq::LongDNA, pos_range::UnitRange{Int})
    for i in reverse(pos_range)
        Base.deleteat!(seq, i)
    end
    return seq
end

# Implemento structure_variation per modifiche maggiori
function structure_variation!(seq::LongDNA{4}, pos::Int, subseq::LongDNA)
    # Test per vedere se è una semplice append
    if pos == lastindex(seq)
        append!(seq, subseq)
        return seq
    # Test per vedere se è una normale insert!
    elseif pos < lastindex(seq)
        insert!(seq, pos, subseq)
        return seq
    # Test per vedere se è una modifica lontana dalla fine della reference
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

# Definiamo i DeltaType
@enum DeltaType DeltaTypeDel DeltaTypeIns DeltaTypeSnp DeltaTypeSV

# Journal Entry
struct JournalEntry
    delta_type::DeltaType    # Tipo della modifica
    position::Int            # Posizione della modifica
    data::Any                # Dati che servono alla modifica ...
end                          # ... (es. una sequenza "GCAT")

# Creo la struttura JournaledString
struct JournaledString
    reference::LongDNA{4}                  # Si può usare anche LongSequence
    journals::Vector{Vector{JournalEntry}} # Vector che definisce le modifiche..
end                                        # ...rispetto alla ref. (i journals)

# Funzione per aggiungere una Delta
function add_delta!(journals, indices::Vector{Int}, 
                    delta_type::DeltaType, position::Int, data::Any)
    for idx in indices
        push!(journals[idx], JournalEntry(delta_type, position, data))
    end
end

function apply_journal(reference::LongDNA{4}, journal::Vector{JournalEntry})
    seq = copy(reference)
    for entry in journal
        # Check sul tipo di modifica
        if entry.delta_type == DeltaTypeDel
            seq = delete_at!(seq, entry.position:(entry.position + 
                   entry.data - 1))  # Il dato dell'entry è la fine del range
        elseif entry.delta_type == DeltaTypeIns
            seq = insert!(seq, entry.position, LongDNA{4}(entry.data))
        elseif entry.delta_type == DeltaTypeSnp
            seq[entry.position] = convert(DNA, entry.data)  
            # Cambio di un singolo nucleotide
        elseif entry.delta_type == DeltaTypeSV
            seq = structure_variation!(seq, entry.position, entry.data)
        end
    end
    return seq
end


# Funzione per stampare le sequenze
function print_sequences(jst::JournaledString)
    for i in 1:length(jst.journals)
        modified_seq = apply_journal(jst.reference, jst.journals[i])
        println("Sequence $i: ", modified_seq)
    end
end


#= Pseudocodice di prova
# Define a reference sequence
reference_seq = LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG")

# Initialize journals for each of the 10 sequences
journals = [Vector{JournalEntry}() for _ in 1:10]

# Create a JournaledString
jst = JournaledString(reference_seq, journals)

# Example: Apply insertions, deletions, SNPs, etc.
println("Prova di: INS e SNP")
add_delta!(jst.journals, [1, 2], DeltaTypeIns, 8, "CGTA")
add_delta!(jst.journals, [4, 9], DeltaTypeSnp, 10, 'C')
add_delta!(jst.journals, [8], DeltaTypeIns, 24, "NNNNN")
print_sequences(jst)

println("\n")
println("Prova di: DELETE (stringa 3 perde 10 caratteri)")
add_delta!(jst.journals, [3], DeltaTypeDel, 1, 10)
print_sequences(jst)

println("\n")
println("Prova di: structure_variation")
add_delta!(jst.journals, [5], DeltaTypeSV, 30, dna"CGTACGTACGTACGTA")
add_delta!(jst.journals, [10], DeltaTypeSV, 18, dna"CGTACGTACGTACGTA")
add_delta!(jst.journals, [7], DeltaTypeSV, 24, dna"NNNNN")
print_sequences(jst)

=#