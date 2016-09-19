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

typealias FourBitAlphs Union{DNANucleotide{4},RNANucleotide{4}}

function count_sites{T<:SiteCase,A<:FourBitAlphs}(::Type{T}, seq::BioSequence{A})
    site_count = 0
    bindex = Seq.bitindex(seq, 1)
    off = Seq.offset(bindex)

    # If the offset of the first element is not zero, then the first integer
    # needs to be shifted / masked.
    if off != 0
        offsetint = seq.data[1] >> off
        site_count += count_sites4(T, offsetint)
    end

    for i in 2:endof(seq.data)
        site_count += count_sites4(T, seq.data[i])
    end

    return site_count
end

#=
function count_sites{T<:SiteCase,A}(::Type{T}, a::BioSequence{A}, b::BioSequence{A})

end
=#
