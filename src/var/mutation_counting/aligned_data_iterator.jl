# Iterating over binary sequence data in an aligned manner
# ========================================================
#
# Not really for export, but for making fast comparrisons between two BioSequences
# which may not have their binary data aligned flush to whole 64 bit integers.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

immutable ShiftedIntsItr{T <: Unsigned}
    vec::Vector{T}
    firstInt::Int
    lastInt::Int
    shift::Int
    revShift::Int

    function ShiftedIntsItr(vec::Vector{T}, firstInt::Int, lastInt::Int, shift::Int)
        @assert firstInt >= 1
        @assert lastInt <= endof(vec)
        ws = sizeof(T) * 8
        @assert shift > 0 && shift <= ws
        return new(vec, firstInt, lastInt, shift, ws - shift)
    end
end

function ShiftedIntsItr(seq::DNASequence)
    bi = Seq.bitindex(seq, 1)
    shift = Seq.offset(bi)
    return ShiftedIntsItr(seq.data,
                          Seq.index(bi),
                          Seq.index(Seq.bitindex(seq, endof(seq))),
                          shift)
end

@inline function Base.start(itr::ShiftedIntsItr)
    return itr.firstInt, itr.firstInt + 1
end

@inline function Base.next(itr::ShiftedIntsItr, state::Tuple{Int64, Int64})
    val = ifelse(state[1] == itr.lastInt,
                 itr.vec[state[1]] >> itr.shift,
                 (itr.vec[state[1]] >> itr.shift) | (itr.vec[state[2]] << itr.revShift))
    return val, (state[2], state[2] + 1)
end

@inline function Base.done(itr::ShiftedIntsItr, state::Tuple{Int64, Int64})
    return state[1] > itr.lastInt
end
