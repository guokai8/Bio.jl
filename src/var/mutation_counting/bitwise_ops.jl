# Bitwise operations for identifying and counting mutations
# =========================================================
#
# Internal bitwise methods used for identifying and counting mutations from the
# 2 and 4 bit encoded representations of DNA and RNA sequences.
#
# These functions probably should not be exported, as it's easy to do things
# very wrong with them if you are not careful. However, they allow the
# counting of mutations between sequences to be wicked fast, whilst taking
# advantage of their compressed representation and avoiding memory expensive
# conversions to large matrices of Nucleotide types, which take a byte of space
# each!
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md


# Nibble operations
# =================

"""
    enumerate_nibbles(abxor::UInt64)

Count the number of set bits, in each nibble of a unsigned 64 bit integer
(With BioSequence 4 bit encoding, each nibble is one RNA or DNA nucleotide).

**This is an internal method and should not be exported.**

E.g. An input of:

0100 0010 0001 0110 1100 1110 1101 1111

Would result in:

0001 0001 0001 0010 0010 0011 0011 0100

This is used to identify different occurances of bit patterns.
"""
@inline function enumerate_nibbles(x::UInt64)
    x = x - ((x >> 1) & 0x5555555555555555)
    return (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
end

"""
    count_zero_nibbles(x::UInt64)

Counts the number of nibbles in a UInt64 `x` that have all their bits unset i.e.
all nibbles of 0000.

**This is an internal method and should not be exported.**

E.g. An input of:

0x0111111111111111

Would give the answer: 1.
"""
@inline function count_zero_nibbles(x::UInt64)
    return 16 - count_ones((x & 0x1111111111111111) |
    (x & 0x2222222222222222) >> 1 |
    (x & 0x4444444444444444) >> 2 |
    (x & 0x8888888888888888) >> 3)
end


# Nibble masking functions
# ------------------------

#= Old version, new version uses fewer operations and is cleaner.
@inline function create_nibble_mask(x::UInt64, value::UInt64)
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

Create a mask for the nibbles (groups of four bits) in a 64 bit integer `x`
that match a given value dictated by the pattern in `value`.

**This is an internal method and should not be exported.**
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

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Gap}, x::UInt64)
    return mask_nibbles(x, 0x0000000000000000)
end

"""
    create_nibble_mask(::Type{Ambiguous}, x::UInt64)

Create a mask of the nibbles in a chunk of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent ambiguous sites.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Ambiguous}, x::UInt64)
    x = enumerate_nibbles(x)
    return mask_nibbles(x, 0x2222222222222222) |
    mask_nibbles(x, 0x3333333333333333) |
    mask_nibbles(x, 0x4444444444444444)
end

"""
    create_nibble_mask(::Type{Pairdel}, x::UInt64)

Create a mask of the nibbles in a chunk of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites that should be
ignored, when counting pairwise mutations between sequences.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Pairdel}, x::UInt64)
    return create_mask(Gap, x) | create_mask(Ambiguous, x)
end

"""
    create_nibble_mask(::Type{Gap}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where either or
both of the two chunks have a gap at a given nibble.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Gap}, a::UInt64, b::UInt64)
    return create_mask(Gap, a) | create_mask(Gap, b)
end

"""
    create_nibble_mask(::Type{Ambiguous}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where either or
both of the two chunks have an ambiguous nucleotide at a given nibble.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Ambiguous}, a::UInt64, b::UInt64)
    return create_mask(Ambiguous, a) | create_mask(Ambiguous, b)
end

"""
    create_nibble_mask(::Type{Pairdel}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where either or
both of the two chunks have a nucleotide symbol that sould be ignored during
pairwise distance computation at a given nibble.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Pairdel}, a::UInt64, b::UInt64)
    return create_mask(Pairdel, a) | create_mask(Pairdel, b)
end


# Bitparallel site counting
# =========================

"""
    count_sites4(::Type{Gap}, x::UInt64)

Count the number of gap sites in a chunk of BioSequence{(DNA|RNA)Nucleotide{4}}
data.

Note that gap sites and empty unused segments of a UInt64 are both 0000, and so
furthur checking of this result would be required in higher level calling
functions.
"""
@inline function count_sites4(::Type{Gap}, x::UInt64)
    return count_zero_nibbles(x)
end

"""
    count_sites4(::Type{Gap}, a::UInt64, b::UInt64)

Count the number of sites in two aligned chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data which contain gap characters.

Note that gap sites and empty unused segments of a UInt64 are both 0000, and so
furthur checking of this result would be required in higher level calling
functions.
"""
@inline function count_sites4(::Type{Gap}, a::UInt64, b::UInt64)
    # Count the gaps in a, count the gaps in b, subtract the number of shared gaps.
    return count_sites4(Gap, a) + count_sites4(Gap, b) - count_sites4(Gap, a | b)
end

"""
    count_sites4(::Type{Ambiguous}, x::UInt64)

Count the number of
ambiguous sites in a chunk of BioSequence{(DNA|RNA)Nucleotide{4}} data.
Ambiuous sites are defined as those with more than one bit set.
Note here gap - 0000 - then is not ambiguous, even though it is a candidate for
pairwise deletion.
"""
@inline function count_sites4(::Type{Ambiguous}, x::UInt64)
    return 16 - count_zero_nibbles(enumerate_nibbles(x) & 0xEEEEEEEEEEEEEEEE)
end

"""
    count_sites4(::Type{Gap}, a::UInt64, b::UInt64)

Count the number of sites in two aligned chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data which contain ambiguous characters.
Ambiuous sites are defined as those with more than one bit set.
Note here gap - 0000 - then is not ambiguous, even though it is a candidate for
pairwise deletion.
"""
@inline function count_sites4(::Type{Ambiguous}, a::UInt64, b::UInt64)
    return 16 - count_zero_nibbles((enumerate_nibbles(a) | enumerate_nibbles(b)) & 0xEEEEEEEEEEEEEEEE)
end






"""
    count_sites4(::Type{Pairdel}, x::UInt64)

An _internal_ function _not for export_, which will count the number of sites in
a chunk of BioSequence{(DNA|RNA)Nucleotide{4}} data that would be ignored in
counts of mutations.
Such sites are defined as those with gaps or ambiguous characters in them.
"""
@inline function count_sites4(::Type{Pairdel}, x::UInt64)
    return count_sites4(Ambiguous, x) + count_sites4(Gap, x)
end

@inline function count_sites4(::Type{Pairdel}, x::UInt64)
    1 + count_sites4(Gap, a, b)
end


"""
    count_sites4(::Type{Match}, a::UInt64, b::UInt64)

An _internal_ function, _not for export_, which will count the number of
mismatches between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data.

**Note:** Ambiguous cases or cases with gaps are ignored and not counted as
matches. For example, 'A' and 'R', or 'A' and '-' will not be counted.
"""
@inline function count_sites4(::Type{Match}, a::UInt64, b::UInt64)
    sharedGaps = count_sites4(Gap, a | b)
    cases = a $ b
    # cases is the result of xoring a and b.
    # if two nucleotides are different, they will contain 1's.
    # if two nucleotides are matches or 2 gaps, they will only have 0's.
    matchGapCount = bitpar_zeros4(cases)
    return matchGapCount - sharedGaps
end

"""
    count_sites4(::Type{Mismatch}, a::UInt64, b::UInt64)

An _internal_ function, _not for export_, which will count the number of
mismatches between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data.

**Note:** Ambiguous cases or cases with gaps are ignored and not counted as
mismatches. For example, 'A' and 'R', or 'A' and '-' will not be counted.
"""
@inline function count_sites4(::Type{Mismatch}, a::UInt64, b::UInt64)
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

@inline function count_sites4(::Type{Ambiguous}, a::UInt64, b::UInt64)
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
    matchBitCount = count_sites4(match)     # 0010 is a normal mismatch.
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
