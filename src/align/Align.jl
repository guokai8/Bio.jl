# Bio.Align
# =========
#
# Module for sequence alignments.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module Align

export
    Operation,
    AlignmentAnchor,
    Alignment,
    AlignedSequence,
    seq2ref,
    ref2seq,
    ismatchop,
    isinsertop,
    isdeleteop,
    cigar,
    OP_MATCH,
    OP_INSERT,
    OP_DELETE,
    OP_SKIP,
    OP_SOFT_CLIP,
    OP_HARD_CLIP,
    OP_PAD,
    OP_SEQ_MATCH,
    OP_SEQ_MISMATCH,
    OP_BACK,
    OP_START,
    # alignment types
    GlobalAlignment,
    SemiGlobalAlignment,
    OverlapAlignment,
    LocalAlignment,
    EditDistance,
    HammingDistance,
    LevenshteinDistance,
    # substitution matrices
    AbstractSubstitutionMatrix,
    SubstitutionMatrix,
    DichotomousSubstitutionMatrix,
    EDNAFULL,
    PAM30,
    PAM70,
    PAM250,
    BLOSUM45,
    BLOSUM50,
    BLOSUM62,
    BLOSUM80,
    BLOSUM90,
    # alignment models
    AbstractScoreModel,
    AffineGapScoreModel,
    AbstractCostModel,
    CostModel,
    # pairwise alignment
    PairwiseAlignment,
    count_matches,
    count_mismatches,
    count_insertions,
    count_deletions,
    count_aligned,
    PairwiseAlignmentResult,
    pairalign,
    score,
    distance,
    alignment,
    hasalignment,
    # HTS
    SAMHeader,
    SAMRecord,
    SAMReader,
    SAMWriter,
    BAMRecord,
    BAMReader,
    BAMWriter,
    BAI,
    AuxDataDict,
    ismapped,
    refname,
    nextrefname,
    refindex,
    nextrefindex,
    leftposition,
    rightposition,
    nextleftposition,
    mappingquality,
    flag,
    templatelength,
    seqname,
    cigar,
    sequence,
    seqlength,
    qualities,
    header,
    isoverlapping,  # re-export from Bio.Intervals
    # SAM flags
    SAM_FLAG_PAIRED,
    SAM_FLAG_PROPER_PAIR,
    SAM_FLAG_UNMAP,
    SAM_FLAG_MUNMAP,
    SAM_FLAG_REVERSE,
    SAM_FLAG_MREVERSE,
    SAM_FLAG_READ1,
    SAM_FLAG_READ2,
    SAM_FLAG_SECONDARY,
    SAM_FLAG_QCFAIL,
    SAM_FLAG_DUP,
    SAM_FLAG_SUPPLEMENTARY

import Bio
using Bio.Seq
using Bio.Intervals
import Bio.StringFields: StringField
import IntervalTrees

include("operations.jl")
include("anchors.jl")
include("alignment.jl")
include("alignedseq.jl")

include("types.jl")
include("submat.jl")
include("models.jl")
include("pairwise/pairalign.jl")
include("hts/hts.jl")

include("deprecated.jl")

end
