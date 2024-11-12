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

# Definisco il current_time
current_time = Ref(0)

# Journal Entry
struct JournalEntry
    delta_type::DeltaType    # Tipo della modifica
    position::Int            # Posizione della modifica
    data::Any                # Dati che servono alla modifica ...
    time::Int                # ... (es. una sequenza "GCAT")
end                          

# Creo la struttura JournaledString
struct JournaledString
    reference::LongDNA{4}                  # Si può usare anche LongSequence
    deltaMap::Vector{Vector{JournalEntry}} # Vector che definisce le modifiche..
end                                        # ...rispetto alla ref. (deltaMap)

# Funzione per aggiungere una Delta
function add_delta!(deltaMap, indices::Vector{Int}, 
                    delta_type::DeltaType, position::Int, data::Any)
    for idx in indices
        push!(deltaMap[idx], JournalEntry(delta_type, position, data,
             current_time[]))
        current_time[] += 1
    end
end

function apply_delta(reference::LongDNA{4}, delta::Vector{JournalEntry})
    seq = copy(reference)
    for entry in delta
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
    for i in 1:length(jst.deltaMap)
        modified_seq = apply_delta(jst.reference, jst.deltaMap[i])
        println("Sequence $i: ", modified_seq)
    end
end

function print_deltas(jst::JournaledString)
    for j in 1:length(jst.deltaMap)
        for i in jst.deltaMap[j]
            println("JournalEntry: $i su stringa indice $j")
        end
    end
end