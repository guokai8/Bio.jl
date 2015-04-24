
# Base interval types and utilities
# ---------------------------------
bitstype 8 Strand

convert(::Type{Strand}, strand::Uint8) = box(Strand, unbox(Uint8, strand))
convert(::Type{Uint8}, strand::Strand) = box(Uint8, unbox(Strand, strand))

const STRAND_NA   = convert(Strand, 0b000)
const STRAND_POS  = convert(Strand, 0b001)
const STRAND_NEG  = convert(Strand, 0b010)
const STRAND_BOTH = convert(Strand, 0b011)

function Base.show(io::IO, strand::Strand)
    if strand == STRAND_NA
        print(io, "(indeterminate strand)")
    elseif strand == STRAND_POS
        print(io, "+")
    elseif strand == STRAND_NEG
        print(io, "-")
    elseif strand == STRAND_BOTH
        print(io, ".")
    else
        print(io, "(undefined strand)")
    end
end


function Base.isless(a::Strand, b::Strand)
    return convert(Uint8, a) < convert(Uint8, b)
end


# Note, just to be clear: this shadows IntervalTrees.Interval
@doc """
A genomic interval specifies interval with some associated metadata.
""" ->
immutable Interval{T} <: IntervalTrees.AbstractInterval{Int64}
    seqname::String
    first::Int64
    last::Int64
    strand::Strand
    metadata::T
end


function IntervalTrees.first(i::Interval)
    return i.first
end


function IntervalTrees.last(i::Interval)
    return i.last
end


function Base.isless{T}(a::Interval{T}, b::Interval{T},
                        seqname_isless::Function=alphanum_isless)
    if a.seqname != b.seqname
        return seqname_isless(a.seqname, b.seqname)
    elseif a.first != b.first
        return a.first < b.first
    elseif a.last != b.last
        return a.last < b.last
    elseif a.strand != b.strand
        return a.strand < b.strand
    else
        return a.metadata < b.metadata
    end
end


function precedes{T}(a::Interval{T}, b::Interval{T},
                     seqname_isless::Function=alphanum_isless)
    return seqname_isless(a.seqname, b.seqname) ||
           a.seqname == b.seqname && a.last < b.first
end


function =={T}(a::Interval{T}, b::Interval{T})
    return a.seqname  == b.seqname &&
           a.first    == b.first &&
           a.last     == b.last &&
           a.strand   == b.strand &&
           a.metadata == b.metadata
end


function isoverlapping{S, T}(a::Interval{S}, b::Interval{T})
    return a.seqname == b.seqname && a.first <= b.last && b.first <= a.last
end


function Base.show(io::IO, i::Interval)
    print(io, i.seqname, ":", i.first, "-", i.last, "    ", i.strand, "    ", i.metadata)
end


@doc """
A type deriving `IntervalStream{T}` must be iterable and produce
Interval{T} objects in sorted order.
""" ->
abstract IntervalStream{T}



@doc """
A comparison function used to sort on numbers within text.

This is useful since sequences are often named things like "chr12" or
"read1234". Treating the numbers as numbers and not text gives a more natural
ordering.

This is similar to the '--version-sort' option in GNU coreutils sort.
""" ->
function alphanum_isless(a::String, b::String)
    i = 1
    j = 1

    # match up to the first digit
    k0 = 0 # position of first digit
    while i <= length(a) && j <= length(b)
        if isdigit(a[i]) && isdigit(b[j])
            k0 = i
            break
        else
            if a[i] != b[j]
                return a[i] < b[j]
            end
        end
        i = nextind(a, i)
        j = nextind(b, j)
    end

    # match numbers
    ka1, kb1 = 0, 0
    while i <= length(a) && isdigit(a[i])
        ka1 = i
        i = nextind(a, i)
    end
    while j <= length(b) && isdigit(b[j])
        kb1 = j
        j = nextind(b, j)
    end

    if ka1 == 0 && kb1 != 0
        return true
    elseif ka1 != 0 && kb1 == 0
        return false
    elseif ka1 != 0 && kb1 != 0
        aval = parse(Int, a[k0:ka1])
        bval = parse(Int, b[k0:kb1])
        if aval != bval
            return aval < bval
        end
    end

    # match suffixes
    while i <= length(a) && j <= length(b)
        if a[i] != b[j]
            return a[i] < b[j]
        end
        i = nextind(a, i)
        j = nextind(b, j)
    end

    return j <= length(b)
end

