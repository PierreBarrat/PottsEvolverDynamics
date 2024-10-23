function filter_X_symbol(aln::AbstractString)
    alphabet = Alphabet(prod(symbols(Alphabet(:aa))) * "X")
    aln = read_fasta(aln; alphabet)
    aln = filter(s -> !in(22, s), aln)
    return aln
end
