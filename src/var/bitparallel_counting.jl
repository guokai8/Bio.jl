abstract SiteCase
immutable Match <: SiteCase end
immutable Mismatch <: SiteCase end
immutable Transition <: SiteCase end
immutable Transversion <: SiteCase end
immutable Gap <: SiteCase end
immutable Ambiguous <: SiteCase end
immutable Pairdel <: SiteCase end



# Internal low-level functions

"""
    bp_count_allzero4()

Recieves a filtered Nucleotide map.

Counts the cases of interest.

A case of interest is anything where all four bits of a nucleotide are not set.

What a case of interest actually means depends on the operations and masks
performed on `x` before this method is called on `x`.

**This is an internal method and should not be exported.**

E.g. An input of:

0001 0000 0001 0010 0010 0011 0011 0100
     ^^^^
Would result in:

0001 0000 0001 0001 0001 0001 0001 0001
     ^^^^
"""
@inline function bp_count_allzero4(x::UInt64)
    return 16 - count_ones((x & 0x1111111111111111) |
    (x & 0x2222222222222222) >> 1 |
    (x & 0x4444444444444444) >> 2 |
    (x & 0x8888888888888888) >> 3)
end

"""
    bp_enumerate4(abxor::UInt64)

Count the number of set bits, in groups of four bits (which represent each abstract Nucleotide).

**This is an internal method and should not be exported.**

E.g. An input of:

0100 0010 0001 0110 1100 1110 1101 1111

Would result in:

0001 0001 0001 0010 0010 0011 0011 0100

This is used to identify different occurances of bit patterns.
"""
@inline function enumerate4(x::UInt64)
    x = x - ((x >> 1) & 0x5555555555555555)
    return (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
end

## Base position masking functions.

#= Old version, new version uses fewer operations and is cleaner.
@inline function mask_nibbles(x::UInt64, value::UInt64)
    x = ~(x $ value)
    x = (x & 0x1111111111111111) &
    ((x >> 1) & 0x1111111111111111) &
    ((x >> 2) & 0x1111111111111111) &
    ((x >> 3) & 0x1111111111111111)
    return x | (x << 1) | (x << 2) | (x << 3)
end
=#

"""
    create_nibble_mask(x::UInt64, value::UInt64)

An internal method, **not for export**, which creates masks for the nibbles
(groups of four bits) in a 64 bit integer `x` that match a given value dictated
by the pattern in `value`.
"""
@inline function create_nibble_mask(value::UInt64, x::UInt64)
    # XOR with the desired values. So matching nibbles will be 0000.
    x $= value
    # Horizontally OR the nibbles.
    x |= (x >> 1)
    x |= (x >> 2)
    # AND removes junk, we then widen x by multiplication and return
    # the inverse.
    x &= 0x1111111111111111
    x *= 15
    return ~x
end


"""
    create_nibble_mask(::Type{Gap}, x::UInt64)

Create a mask of the nibbles in a chunk of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent gaps.
"""
@inline function create_nibble_mask(::Type{Gap}, x::UInt64)
    return mask_nibbles(x, 0x0000000000000000)
end

"""
    create_nibble_mask(::Type{Ambiguous}, x::UInt64)

Create a mask of the nibbles in a chunk of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent ambiguous sites.
"""
@inline function create_nibble_mask(::Type{Ambiguous}, x::UInt64)
    c = enumerate4(x)
    return mask_nibbles(c, 0x2222222222222222) |
    mask_nibbles(c, 0x3333333333333333) |
    mask_nibbles(c, 0x4444444444444444)
end

"""
    create_nibble_mask(::Type{Pairdel}, x::UInt64)

Create a mask of the nibbles in a chunk of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites that should be
ignored, when counting pairwise mutations between sequences.
"""
@inline function create_nibble_mask(::Type{Pairdel}, x::UInt64)
    return create_mask(Gap, x) | create_mask(Ambiguous, x)
end

@inline function create_nibble_mask(::Type{Gap}, a::UInt64, b::UInt64)
    return create_mask(Gap, a) | create_mask(Gap, b)
end

@inline function create_nibble_mask(::Type{Ambiguous}, a::UInt64, b::UInt64)
    return create_mask(Ambiguous, a) | create_mask(Ambiguous, b)
end

@inline function create_nibble_mask(::Type{Pairdel}, a::UInt64, b::UInt64)
    return create_mask(Pairdel, a) | create_mask(Pairdel, b)
end















"""
    bitpar_count4(::Type{Gap}, x::UInt64)

An _internal_ function _not for export_, which will count the number of
gap sites in a chunk of BioSequence{(DNA|RNA)Nucleotide{4}} data.
Note that gap sites and empty unused segments of a UInt64 are both 0000,
and so furthur checking of this result would be required in higher level calling
functions.
"""
@inline function bitpar_count4(::Type{Gap}, x::UInt64)
    return bp_count_allzero4(x)
end

@inline function _mask(::Type{Gap}, x::UInt64)
    return bp_enumerate4(~x)
end

@inline function _mask(::Type{Gap}, a::UInt64, b::UInt64)
    return _mask(Gap, a) | _mask(Gap, b)
end





"""
    bitpar_count4(::Type{Gap}, a::UInt64, b::UInt64)

An _internal_ function _not for export_, which will count the number of
**positions** in two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data which
contain gap characters.
"""
@inline function bitpar_count4(::Type{Gap}, a::UInt64, b::UInt64)
    # Count the gaps in a, count the gaps in b, subtract the number of shared gaps.
    return bitpar_count4(Gap, a) + bitpar_count4(Gap, b) - bitpar_count4(Gap, a | b)
end

"""
    bitpar_count4(::Type{Ambiguous}, x::UInt64)

An _internal_ function _not for export_, which will count the number of
ambiguous sites in a chunk of BioSequence{(DNA|RNA)Nucleotide{4}} data.
Ambiuous sites are defined as those with more than one bit set.
Note here gap - 0000 - then is not ambiguous, even though it is a candidate for
pairwise deletion.
"""
@inline function bitpar_count4(::Type{Ambiguous}, x::UInt64)
    return 16 - bp_count_allzero4(bp_enumerate4(x) & 0xEEEEEEEEEEEEEEEE)
end

@inline function bitpar_count4(::Type{Ambiguous}, a::UInt64, b::UInt64)

end

"""
    bitpar_count4(::Type{Pairdel}, x::UInt64)

An _internal_ function _not for export_, which will count the number of sites in
a chunk of BioSequence{(DNA|RNA)Nucleotide{4}} data that would be ignored in
counts of mutations.
Such sites are defined as those with gaps or ambiguous characters in them.
"""
@inline function bitpar_count4(::Type{Pairdel}, x::UInt64)
    return bitpar_count4(Ambiguous, x) + bitpar_count4(Gap, x)
end

@inline function bitpar_count4(::Type{Pairdel}, x::UInt64)
    + bitpar_count4(Gap, a, b)
end


"""
    bitpar_count4(::Type{Match}, a::UInt64, b::UInt64)

An _internal_ function, _not for export_, which will count the number of
mismatches between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data.

**Note:** Ambiguous cases or cases with gaps are ignored and not counted as
matches. For example, 'A' and 'R', or 'A' and '-' will not be counted.
"""
@inline function bitpar_count4(::Type{Match}, a::UInt64, b::UInt64)
    sharedGaps = bitpar_count4(Gap, a | b)
    cases = a $ b
    # cases is the result of xoring a and b.
    # if two nucleotides are different, they will contain 1's.
    # if two nucleotides are matches or 2 gaps, they will only have 0's.
    matchGapCount = bitpar_zeros4(cases)
    return matchGapCount - sharedGaps
end

"""
    bitpar_count4(::Type{Mismatch}, a::UInt64, b::UInt64)

An _internal_ function, _not for export_, which will count the number of
mismatches between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data.

**Note:** Ambiguous cases or cases with gaps are ignored and not counted as
mismatches. For example, 'A' and 'R', or 'A' and '-' will not be counted.
"""
@inline function bitpar_count4(::Type{Mismatch}, a::UInt64, b::UInt64)
    cases = a $ b
    # cases is the result of xoring a and b.
    # If two nucleotides are different, they will contain 1's.
    # If two nucleotides are matches or 2 gaps, they will only have 0's.
    # Unambiguous mismatches always have 2 set bits, and we can explot this.
    enumeratedCases = bitpar_enumerate4(cases)
    # enumeratedCases contains the number of set bits, for each position.
    # When a position is 0000, it either represents a match or a gap.
    # A normal mismatch is 0010, anything ambiguous will not be 0010.
    matchMismatchGapCount = bitpar_zeros4(enumeratedCases & 0xDDDDDDDDDDDDDDDD)
    # The enumeratedCases are filtered by masking with 0xDDDDDD...
    # this results in clear mismatches i.e. 0010 being masked to 0000.
    # These 0000 cases are then counted by bitpar_zeros4 to get the number of
    # cases that are clear matches, clear mismatches, and those that are gaps.
    # To get the number of mismatches, we now have to enumerate the number of
    # matches or gaps, and subtract that from matchMismatchGapCount.
    matchGapCount = bitpar_zeros4(cases)
    return matchMismatchGapCount - matchGapCount
end

@inline function bitpar_count4(::Type{Ambiguous}, a::UInt64, b::UInt64)
    cases = a $ b
    # cases is the result of xoring a and b.
    # If two nucleotides are different, they will contain 1's.
    # If two nucleotides are matches or 2 gaps, they will only have 0's.
    # Unambiguous mismatches always have 2 set bits, and we can explot this.
    enumeratedCases = bitpar_enumerate4(cases)
    # enumeratedCases contains the number of set bits, for each position.
    # When a position is 0000, it either represents a match or a gap.
    # A normal mismatch is 0010, anything ambiguous will not be 0010.
    matchAmbiguousGapCount = bitpar_zeros4(enumeratedCases & 0x2222222222222222)
    # The enumeratedCases are filtered by masking with 0x2222222...
    # this results in ambiguous cases i.e. not 0010, being masked to 0000.
    # These 0000 cases are then counted by bitpar_zeros4 to get the number of
    # cases that are clear matches, ambiguous cases, and those that are gaps.
    # To get the number of ambiguous cases, we now have to enumerate the number
    # of matches or gaps, and subtract that from matchAmbiguousGapCount.
    matchGapCount = bitpar_zeros4(cases)
    return matchAmbiguousGapCount - matchGapCount
end








function bitpar_mismatches_gaps(a::UInt64, b::UInt64)
    gapCount = gapcount(a, b)
    match = a $ b
    matchGapCount = bitpar_zeros4(match)
    # matchBitCount contains the number of set bits, for each position.
    # when a position is 0000, it either represents a match or a gap.
    # A normal mismatch contains only two bits, now we will count them.
    matchBitCount = bitpar_count4(match)     # 0010 is a normal mismatch.
    filteredMatches = bitpar_zeros4(matchBitCount & 0xDDDDDDDDDDDDDDDD)
    # filteredMatches represents the cases where the result was not 0010.
    return filteredMatches - matchGapCount
end

@inline function transition(a::UInt64, b::UInt64)
    gapCount = gapcount(a, b)
    abxor = a $ b
    ctTransitionFilter = abxor $ 0xAAAAAAAAAAAAAAAA
    gaTransitionFilter = abxor $ 0x5555555555555555
    ctTransitionCount = bitpar_zeros4(ctTransitionFilter)
    gaTransitionCount = bitpar_zeros4(gaTransitionFilter)
    gapsOrMatchesCount = bitpar_zeros4(abxor)
    return ctTransitionCount + gaTransitionCount
end
