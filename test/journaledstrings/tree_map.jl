using DataStructures
include("journaledstrings.jl")

Jtree = AVLTree{Int64}()

Jmap = DefaultDict{Int64, JournalEntry}()
