
# Iterator to iterate over sequences
immutable AlignedDataItr
    seq::DNASequence
    frontMask::UInt64
    backMask::UInt64
    firstInt::Int
    lastInt::Int
    off::Int
    revoff::Int
end

function aligned_data(seq::DNASequence)
    bi = Seq.bitindex(seq, 1)
    off = Seq.offset(bi)
    bm = Seq.mask(off)
    fm = ~bm
    return AlignedDataItr(seq, fm, bm,
                          Seq.index(bi),
                          Seq.index(Seq.bitindex(seq, endof(seq))),
                          off,
                          64 - off)
end

@inline Base.start(itr::AlignedDataItr) = (firstInt, firstInt + 1)

@inline function Base.next(itr::AlignedDataItr, state::Tuple{Int64, Int64})
    # There may not even be any need of frontMask and backMask.
    # val = ((itr.seq.data[1] & itr.frontMask) >> itr.off) |
    #      ((itr.seq.data[2] & itr.backMask) << (64 - itr.revoff))

    val = (itr.seq.data[state[1]] >> itr.off)
    val |= (itr.seq.data[state[2]] << itr.revoff)

    return val, state
end

@inline function Base.done(itritr::AlignedDataItr, state::Tuple{Int64, Int64})

end



immutable AlignedIteratorState
    at::UInt64
    first::UInt64
    last::UInt64
    len::UInt64
    firstoffset::UInt8
    lastoffset::UInt8
end

function Base.start(iter::AlignedIterator)
    return AlignedIteratorState(1,
    Seq.index(Seq.bitindex(iter.seq, 1)),
    Seq.index(Seq.bitindex(iter.seq, endof(iter.seq))),
    (Seq.index(Seq.bitindex(iter.seq, endof(iter.seq))) - Seq.index(Seq.bitindex(iter.seq, 1))) + 1,
    Seq.offset(Seq.bitindex(iter.seq, 1)),
    Seq.offset(Seq.bitindex(iter.seq, endof(iter.seq)))
    )
end

function next(iter::AlignedIterator, state)
    if state.firstoffset == 0 && state.lastoffset == 60
        return (iter.seq.data[state.first + state.at], AlignedIteratorState(
        state.at + 1,
        state.first,
        state.last,
        state.len,
        state.firstoffset,
        state.lastoffset))
    elseif state.firstoffset != 0 && state.lastoffset == 60

    elseif state.firstoffset == 0 && state.lastoffset != 60

    elseif state.firstoffset != 0 && state.lastoffset != 60

    end
end

function done(iter, state)
    return state.at > state.len
end
