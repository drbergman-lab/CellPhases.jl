# The PhysiCellOutput extension, on PhysiCell-format output written at test time: 150 frames, 10
# minutes apart, of rgc cells carrying a plane wave in `notch`, plus fixed "layer" cells. Cell 1
# turns into a layer cell for frames 70–71 and cell 2 is dead in frame 100, which splits their
# tracks. The PhysiCellModelManager extension only resolves simulations to output folders and is
# not loaded here.

import PhysiCellOutput
include("physicell_output.jl")

@testset "PhysiCellOutput extension" begin
    rgc = agentField(; n = 300, L = 250.0, nframes = 150, step = 0.5, rng = Xoshiro(21))
    notch = planeWave(rgc.positions; dt = 10, period = 180, speed = 0.5, direction = 90, rng = Xoshiro(22))
    dir = writePhysiCellOutput(mktempdir(), (; rgc.positions, values = notch, rgc.ids);
                               layer = Dict(70 => [1], 71 => [1]), dead = Dict(100 => [2]))
    ext = Base.get_extension(CellPhases, :CellPhasesPhysiCellOutputExt)
    sequence = PhysiCellOutput.PhysiCellSequence(dir)

    f = ext._frames(sequence; value = :notch, cell_types = "rgc")
    @test f.dt == 10.0
    @test length(f.ids) == 150
    @test f.ids[1] == 1:300
    @test 1 ∉ f.ids[70] && 2 ∉ f.ids[100] && 2 ∈ f.ids[99]
    @test f.values[5] ≈ notch[5]
    @test [x[2] for x in f.positions[5]] ≈ [x[2] for x in rgc.positions[5]]
    @test length(first(ext._frames(sequence; value = :notch, dims = 3).positions[1])) == 3
    both = ext._frames(sequence; value = "notch", cell_types = ["rgc", "layer"])
    @test length(both.ids[1]) == 330
    @test length(ext._frames(sequence; value = :notch, include_dead = true).ids[100]) == 330

    # end to end: the same phases as the array method, and the wave recovered
    p = cellPhases(sequence; value = :notch, cell_types = "rgc", band = (100, 300), minlength = 600)
    ref = cellPhases(f.positions, f.values, f.ids; dt = 10, band = (100, 300), minlength = 600)
    @test isequal(p.phase, ref.phase)
    @test [tr.frames for tr in p.tracks if tr.id == 1] == [19:51, 90:132]
    @test [tr.frames for tr in p.tracks if tr.id == 2] == [19:81]     # frames 101–150 too short
    w = PhaseWaves(p; r = 30, frames = 1:4:150)
    good = w.r2 .>= 0.5
    @test rad2deg(angle(mean(cis.(w.direction[good])))) ≈ 90 atol = 3
    @test median(filter(!isnan, w.speed[good])) ≈ 0.5 rtol = 0.05

    @test_throws ArgumentError cellPhases(sequence; value = :notch, cell_types = "rcg", band = (100, 300))
    @test_throws ArgumentError cellPhases(sequence; value = :notch, dt = 10, band = (100, 300))
    @test_throws ArgumentError ext._frames(sequence; value = :notch, dims = 4)
    @test_throws Exception ext._frames(sequence; value = :not_a_label)
end
