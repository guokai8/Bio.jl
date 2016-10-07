# Types and methods for counting mutations
# ========================================
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

abstract SiteCase
immutable Conserved <: SiteCase end
immutable Mutated <: SiteCase end
immutable Transition <: SiteCase end
immutable Transversion <: SiteCase end
immutable Gap <: SiteCase end
immutable Ambiguous <: SiteCase end
immutable Pairdel <: SiteCase end

include("bitwise_ops.jl")
include("aligned_data_iterator.jl")

typealias FourBitAlphs Union{DNAAlphabet{4},RNAAlphabet{4}}

function count_sites{A<:FourBitAlphs}(::Type{Ambiguous}, seq::BioSequence{A})
    site_count = 0
    bindex = Seq.bitindex(seq, 1)
    firstoff = Seq.offset(bindex)

    # If the offset of the first element is not zero, then the first integer
    # needs to be shifted / masked.
    if firstoff != 0
        offsetint = seq.data[1] >> firstoff
        site_count += count_sites4(T, offsetint)
        for i in 2:endof(seq.data)
            site_count += count_sites4(Ambiguous, seq.data[i])
        end
    else
        for i in 1:endof(seq.data)
            site_count += count_sites4(Ambiguous, seq.data[i])
        end
    end
    return site_count
end








function count_sites{T<:SiteCase,A}(::Type{T}, a::BioSequence{A}, b::BioSequence{A})

    @assert length(a) == length(b)


end
