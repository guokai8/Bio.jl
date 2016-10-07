
# Iterator to iterate over binary sequence data in aligned manner.
immutable AlignedDataItr
    seq::DNASequence
    lastInt::Int
    off::Int
    revoff::Int
end

function aligned_data(seq::DNASequence)
    off = Seq.offset(Seq.bitindex(seq, 1))
    return AlignedDataItr(seq, Seq.index(Seq.bitindex(seq, endof(seq))), off, 64 - off)
end

@inline function Base.start(itr::AlignedDataItr)
    i = Seq.index(Seq.bitindex(itr.seq, 1))
    return (i, i + 1)
end

@inline function Base.next(itr::AlignedDataItr, state::Tuple{Int64, Int64})
    val = itr.seq.data[state[1]] >> itr.off
    val |= ifelse(state[1] == itr.lastInt,
                  UInt64(0),
                  itr.seq.data[state[2]] << itr.revoff)
    return val, (state[2], state[2] + 1)
end

@inline function Base.done(itr::AlignedDataItr, state::Tuple{Int64, Int64})
    return state[1] > itr.lastInt
end
