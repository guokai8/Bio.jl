module TestVar

using Base.Test

using Bio.Seq
using Bio.Var

@testset "Var" begin

    @testset "Site counting and identification" begin

        @testset "Bit parallel internals" begin

            # 4 consecutive Conservedes, 4 consecutive misConservedes, 4 consecutive ambigs,
            # 4 consecutive gaps.
            aseq = dna"ATCGATCGARMGA-M-"
            bseq = dna"ATCGTCGARMWY-T--"
            cseq = dna"ATCGMRWSYKVHDBN-"
            a = aseq.data[1]
            b = bseq.data[1]
            c = cseq.data[1]

            @testset "Counting zeros" begin
                @test Var.count_zero_nibbles(0x0000000000000000) == 16
                @test Var.count_zero_nibbles(0xF004020000403010) == 10
            end

            @testset "Enumerating nibbles" begin
                @test Var.enumerate_nibbles(0x0000000000000000) == 0x0000000000000000
                @test Var.enumerate_nibbles(0xF004020000403010) == 0x4001010000102010
            end

            @testset "Masking nibbles" begin
                @test Var.create_nibble_mask(Gap, a) == 0xF0F0000000000000
                @test Var.create_nibble_mask(Gap, b) == 0xFF0F000000000000
                @test Var.create_nibble_mask(Gap, a, b) == 0xFFFF000000000000
                @test Var.create_nibble_mask(Ambiguous, a) == 0x0F000FF000000000
                @test Var.create_nibble_mask(Ambiguous, b) == 0x0000FFFF00000000
                @test Var.create_nibble_mask(Ambiguous, a, b) == 0x0F00FFFF00000000
                @test Var.create_nibble_mask(Pairdel, a) == 0xFFF00FF000000000
                @test Var.create_nibble_mask(Pairdel, b) == 0xFF0FFFFF00000000
                @test Var.create_nibble_mask(Pairdel, a, b) == 0xFFFFFFFF00000000
                @test Var.create_nibble_mask(Conserved, a, b) == 0x000000000000FFFF
                @test Var.create_nibble_mask(Mutated, a, b) == 0x00000000FFFF0000
            end

            @testset "Case counting" begin

                @testset "Gaps" begin
                    @test Var.count_sites4(Gap, a) == 2
                    @test Var.count_sites4(Gap, b) == 3
                    @test Var.count_sites4(Gap, a | b) == Var.count_sites4(Gap, b | a) == 1
                    @test Var.count_sites4(Gap, a, b) == Var.count_sites4(Gap, b, a) == 4
                    @test Var.count_sites4(Gap, a | c) == Var.count_sites4(Gap, c | a) == 1
                    @test Var.count_sites4(Gap, a, c) == Var.count_sites4(Gap, c, a) == 2
                    @test Var.count_sites4(Gap, b, c) == Var.count_sites4(Gap, c, b) == 3
                    @test Var.count_sites4(Gap, b | c) == Var.count_sites4(Gap, c | b) == 1
                end

                @testset "Ambiguities" begin
                    @test Var.count_sites4(Ambiguous, c) == 11
                    @test Var.count_sites4(Ambiguous, a) == 3
                    @test Var.count_sites4(Ambiguous, b) == 4
                    @test Var.count_sites4(Ambiguous, a, b) == Var.count_sites4(Ambiguous, b, a) == 5
                    @test Var.count_sites4(Ambiguous, a, c) == Var.count_sites4(Ambiguous, c, a) == 11
                    @test Var.count_sites4(Ambiguous, b, c) == Var.count_sites4(Ambiguous, c, b) == 11
                end

                @testset "Pairdel" begin
                    @test Var.count_sites4(Pairdel, a) == 5
                    @test Var.count_sites4(Pairdel, b) == 7
                    @test Var.count_sites4(Pairdel, c) == 12
                    @test Var.count_sites4(Pairdel, a, b) == Var.count_sites4(Pairdel, b, a) == 8
                    @test Var.count_sites4(Pairdel, a, c) == Var.count_sites4(Pairdel, c, a) == 12
                    @test Var.count_sites4(Pairdel, b, c) == Var.count_sites4(Pairdel, c, b) == 12
                end

                @testset "Conserved" begin
                    @test Var.count_sites4(Conserved, a, b) == Var.count_sites4(Conserved, b, a) == 4
                    @test Var.count_sites4(Conserved, a, c) == Var.count_sites4(Conserved, c, a) == 4
                    @test Var.count_sites4(Conserved, b, c) == Var.count_sites4(Conserved, c, b) == 4
                end

                @testset "Mutated" begin
                    @test Var.count_sites4(Mutated, a, b) == Var.count_sites4(Mutated, b, a) == 4
                end
            end
        end
    end

    @testset "Old functionality" begin
        @testset "Counting mutations" begin

            # Create a 20bp test DNA sequence pair containing every possible transition (4),
            # every possible transversion (8), and 2 gapped sites and 2 ambiguous sites.
            # This leaves 4 sites non-mutated/conserved.
            dnas = [dna"ATTG-ACCTGGNTTTCCGAA", dna"A-ACAGAGTATACRGTCGTC"]
            m1 = seqmatrix(dnas, :seq)

            rnas = [rna"AUUG-ACCUGGNUUUCCGAA", rna"A-ACAGAGUAUACRGUCGUC"]
            m2 = seqmatrix(rnas, :seq)

            @test count_mutations(AnyMutation, dnas) == count_mutations(AnyMutation, rnas) == ([12], [16])
            @test count_mutations(AnyMutation, m1) == count_mutations(AnyMutation, m2) == ([12], [16])
            @test count_mutations(TransitionMutation, dnas) == count_mutations(TransitionMutation, rnas) == ([4], [16])
            @test count_mutations(TransitionMutation, m1) == count_mutations(TransitionMutation, m2) == ([4], [16])
            @test count_mutations(TransversionMutation, dnas) == count_mutations(TransversionMutation, rnas) == ([8], [16])
            @test count_mutations(TransversionMutation, m1) == count_mutations(TransversionMutation, m2) == ([8], [16])
            @test count_mutations(TransitionMutation, TransversionMutation, dnas) == count_mutations(TransitionMutation, TransversionMutation, rnas) == ([4], [8], [16])
            @test count_mutations(TransitionMutation, TransversionMutation, m1) == count_mutations(TransitionMutation, TransversionMutation, m2) == ([4], [8], [16])
            @test count_mutations(TransversionMutation, TransitionMutation, dnas) == count_mutations(TransversionMutation, TransitionMutation, rnas) == ([4], [8], [16])
            @test count_mutations(TransversionMutation, TransitionMutation, m1) == count_mutations(TransversionMutation, TransitionMutation, m2) == ([4], [8], [16])

            ans = Bool[false, false, true, true, false, true, true, true, false, true, true, false, true, false, true, true, false, false, true, true]
            @test flagmutations(AnyMutation, m1)[1][:,1] == ans
            @test flagmutations(AnyMutation, m2)[1][:,1] == ans
        end

        @testset "Distance Computation" begin

            dnas1 = [dna"ATTG-ACCTGGNTTTCCGAA", dna"A-ACAGAGTATACRGTCGTC"]
            m1 = seqmatrix(dnas1, :seq)

            dnas2 = [dna"attgaacctggntttccgaa", dna"atacagagtatacrgtcgtc"]
            m2 = seqmatrix(dnas2, :seq)

            @test distance(Count{AnyMutation}, dnas1) == ([12], [16])
            @test distance(Count{TransitionMutation}, dnas1) == ([4], [16])
            @test distance(Count{TransversionMutation}, dnas1) == ([8], [16])
            @test distance(Count{Kimura80}, dnas1) == ([4], [8], [16])
            @test distance(Count{AnyMutation}, m1) == ([12], [16])
            @test distance(Count{TransitionMutation}, m1) == ([4], [16])
            @test distance(Count{TransversionMutation}, m1) == ([8], [16])
            @test distance(Count{Kimura80}, m1) == ([4], [8], [16])

            @test distance(Count{AnyMutation}, dnas2) == ([12], [18])
            @test distance(Count{TransitionMutation}, dnas2) == ([4], [18])
            @test distance(Count{TransversionMutation}, dnas2) == ([8], [18])
            @test distance(Count{Kimura80}, dnas2) == ([4], [8], [18])
            @test distance(Count{AnyMutation}, m2) == ([12], [18])
            @test distance(Count{TransitionMutation}, m2) == ([4], [18])
            @test distance(Count{TransversionMutation}, m2) == ([8], [18])
            @test distance(Count{Kimura80}, m2) == ([4], [8], [18])

            @test distance(Proportion{AnyMutation}, dnas1) == ([(12 / 16)], [16])
            @test distance(Proportion{TransitionMutation}, dnas1) == ([(4 / 16)], [16])
            @test distance(Proportion{TransversionMutation}, dnas1) == ([(8 / 16)], [16])
            @test distance(Proportion{AnyMutation}, m1) == ([(12 / 16)], [16])
            @test distance(Proportion{TransitionMutation}, m1) == ([(4 / 16)], [16])
            @test distance(Proportion{TransversionMutation}, m1) == ([(8 / 16)], [16])

            @test distance(Proportion{AnyMutation}, dnas2) == ([(12 / 18)], [18])
            @test distance(Proportion{TransitionMutation}, dnas2) == ([(4 / 18)], [18])
            @test distance(Proportion{TransversionMutation}, dnas2) == ([(8 / 18)], [18])
            @test distance(Proportion{AnyMutation}, m2) == ([(12 / 18)], [18])
            @test distance(Proportion{TransitionMutation}, m2) == ([(4 / 18)], [18])
            @test distance(Proportion{TransversionMutation}, m2) == ([(8 / 18)], [18])

            @test distance(JukesCantor69, dnas1) == ([Inf], [Inf]) # Returns infinity as 12/16 is 0.75 - mutation saturation.
            @test distance(JukesCantor69, m1) == ([Inf], [Inf])

            @test round(distance(JukesCantor69, dnas2)[1][1], 3) == 1.648
            @test round(distance(JukesCantor69, dnas2)[2][1], 3) == 1
            @test round(distance(JukesCantor69, m2)[1][1], 3) == 1.648
            @test round(distance(JukesCantor69, m2)[2][1], 3) == 1

            @test round(distance(Kimura80, dnas2)[1][1], 3) == 1.648
            @test round(distance(Kimura80, dnas2)[2][1], 3) == 1
            @test round(distance(Kimura80, m2)[1][1], 3) == 1.648
            @test round(distance(Kimura80, m2)[2][1], 3) == 1
        end
    end
end

@testset "Counting mutations" begin

    # Create a 20bp test DNA sequence pair containing every possible transition (4),
    # every possible transversion (8), and 2 gapped sites and 2 ambiguous sites.
    # This leaves 4 sites non-mutated/conserved.
    dnas = [dna"ATTG-ACCTGGNTTTCCGAA", dna"A-ACAGAGTATACRGTCGTC"]
    m1 = seqmatrix(dnas, :seq)

    rnas = [rna"AUUG-ACCUGGNUUUCCGAA", rna"A-ACAGAGUAUACRGUCGUC"]
    m2 = seqmatrix(rnas, :seq)

    @test count_mutations(AnyMutation, dnas) == count_mutations(AnyMutation, rnas) == ([12], [16])
    @test count_mutations(AnyMutation, m1) == count_mutations(AnyMutation, m2) == ([12], [16])
    @test count_mutations(TransitionMutation, dnas) == count_mutations(TransitionMutation, rnas) == ([4], [16])
    @test count_mutations(TransitionMutation, m1) == count_mutations(TransitionMutation, m2) == ([4], [16])
    @test count_mutations(TransversionMutation, dnas) == count_mutations(TransversionMutation, rnas) == ([8], [16])
    @test count_mutations(TransversionMutation, m1) == count_mutations(TransversionMutation, m2) == ([8], [16])
    @test count_mutations(TransitionMutation, TransversionMutation, dnas) == count_mutations(TransitionMutation, TransversionMutation, rnas) == ([4], [8], [16])
    @test count_mutations(TransitionMutation, TransversionMutation, m1) == count_mutations(TransitionMutation, TransversionMutation, m2) == ([4], [8], [16])
    @test count_mutations(TransversionMutation, TransitionMutation, dnas) == count_mutations(TransversionMutation, TransitionMutation, rnas) == ([4], [8], [16])
    @test count_mutations(TransversionMutation, TransitionMutation, m1) == count_mutations(TransversionMutation, TransitionMutation, m2) == ([4], [8], [16])

    ans = Bool[false, false, true, true, false, true, true, true, false, true, true, false, true, false, true, true, false, false, true, true]
    @test flagmutations(AnyMutation, m1)[1][:,1] == ans
    @test flagmutations(AnyMutation, m2)[1][:,1] == ans


end

@testset "Distance Computation" begin

    dnas1 = [dna"ATTG-ACCTGGNTTTCCGAA", dna"A-ACAGAGTATACRGTCGTC"]
    m1 = seqmatrix(dnas1, :seq)

    dnas2 = [dna"attgaacctggntttccgaa",
             dna"atacagagtatacrgtcgtc"]
    dnas3 = [dna"attgaacctgtntttccgaa",
             dna"atagaacgtatatrgccgtc"]
    m2 = seqmatrix(dnas2, :seq)

    @test distance(Count{AnyMutation}, dnas1) == ([12], [16])
    @test distance(Count{TransitionMutation}, dnas1) == ([4], [16])
    @test distance(Count{TransversionMutation}, dnas1) == ([8], [16])
    @test distance(Count{Kimura80}, dnas1) == ([4], [8], [16])
    @test distance(Count{AnyMutation}, m1) == ([12], [16])
    @test distance(Count{TransitionMutation}, m1) == ([4], [16])
    @test distance(Count{TransversionMutation}, m1) == ([8], [16])
    @test distance(Count{Kimura80}, m1) == ([4], [8], [16])

    @test distance(Count{AnyMutation}, dnas2, 5, 5)[1][:] == [2, 4, 3, 3]
    @test distance(Count{AnyMutation}, dnas2, 5, 5)[2][:] == [5, 5, 3, 5]
    @test distance(Count{TransitionMutation}, dnas2, 5, 5)[1][:] == [0, 2, 1, 1]
    @test distance(Count{TransitionMutation}, dnas2, 5, 5)[2][:] == [5, 5, 3, 5]
    @test distance(Count{TransversionMutation}, dnas2, 5, 5)[1][:] == [2, 2, 2, 2]
    @test distance(Count{TransversionMutation}, dnas2, 5, 5)[2][:] == [5, 5, 3, 5]
    @test distance(Count{Kimura80}, dnas1) == ([4], [8], [16])

    @test distance(Count{AnyMutation}, dnas2) == ([12], [18])
    @test distance(Count{TransitionMutation}, dnas2) == ([4], [18])
    @test distance(Count{TransversionMutation}, dnas2) == ([8], [18])
    @test distance(Count{Kimura80}, dnas2) == ([4], [8], [18])
    @test distance(Count{AnyMutation}, m2) == ([12], [18])
    @test distance(Count{TransitionMutation}, m2) == ([4], [18])
    @test distance(Count{TransversionMutation}, m2) == ([8], [18])
    @test distance(Count{Kimura80}, m2) == ([4], [8], [18])

    d = distance(Proportion{AnyMutation}, dnas2, 5, 5)
    a = [0.4, 0.8, 1.0, 0.6]
    for i in 1:length(d[1])
        @test_approx_eq_eps d[1][i] a[i] 1e-4
    end
    @test d[2][:] == [5, 5, 3, 5]
    d = distance(Proportion{TransitionMutation}, dnas2, 5, 5)
    a = [0.0, 0.4, 0.333333, 0.2]
    for i in 1:length(d[1])
        @test_approx_eq_eps d[1][i] a[i] 1e-4
    end
    @test d[2][:] == [5, 5, 3, 5]
    d = distance(Proportion{TransversionMutation}, dnas2, 5, 5)
    a = [0.4, 0.4, 0.666667, 0.4]
    for i in 1:length(d[1])
        @test_approx_eq_eps d[1][i] a[i] 1e-4
    end
    @test d[2][:] == [5, 5, 3, 5]

    @test distance(Proportion{AnyMutation}, dnas1) == ([(12 / 16)], [16])
    @test distance(Proportion{TransitionMutation}, dnas1) == ([(4 / 16)], [16])
    @test distance(Proportion{TransversionMutation}, dnas1) == ([(8 / 16)], [16])
    @test distance(Proportion{AnyMutation}, m1) == ([(12 / 16)], [16])
    @test distance(Proportion{TransitionMutation}, m1) == ([(4 / 16)], [16])
    @test distance(Proportion{TransversionMutation}, m1) == ([(8 / 16)], [16])

    @test distance(Proportion{AnyMutation}, dnas2) == ([(12 / 18)], [18])
    @test distance(Proportion{TransitionMutation}, dnas2) == ([(4 / 18)], [18])
    @test distance(Proportion{TransversionMutation}, dnas2) == ([(8 / 18)], [18])
    @test distance(Proportion{AnyMutation}, m2) == ([(12 / 18)], [18])
    @test distance(Proportion{TransitionMutation}, m2) == ([(4 / 18)], [18])
    @test distance(Proportion{TransversionMutation}, m2) == ([(8 / 18)], [18])

    @test distance(JukesCantor69, dnas1) == ([Inf], [Inf]) # Returns infinity as 12/16 is 0.75 - mutation saturation.
    @test distance(JukesCantor69, m1) == ([Inf], [Inf])

    @test round(distance(JukesCantor69, dnas2)[1][1], 3) == 1.648
    @test round(distance(JukesCantor69, dnas2)[2][1], 3) == 1
    @test round(distance(JukesCantor69, m2)[1][1], 3) == 1.648
    @test round(distance(JukesCantor69, m2)[2][1], 3) == 1
    @test_throws DomainError distance(JukesCantor69, dnas2, 5, 5)
    d = distance(JukesCantor69, dnas3, 5, 5)
    a = [0.232616, 0.571605, 0.44084, 0.571605]
    v = [0.0595041, 0.220408, 0.24, 0.220408]
    for i in 1:length(d[1])
        @test_approx_eq_eps d[1][i] a[i] 1e-5
        @test_approx_eq_eps d[2][i] v[i] 1e-5
    end

    @test round(distance(Kimura80, dnas2)[1][1], 3) == 1.648
    @test round(distance(Kimura80, dnas2)[2][1], 3) == 1
    @test round(distance(Kimura80, m2)[1][1], 3) == 1.648
    @test round(distance(Kimura80, m2)[2][1], 3) == 1

end

end # module TestVar
