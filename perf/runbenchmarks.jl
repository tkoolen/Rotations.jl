using Rotations
using BenchmarkTools
import Base.Iterators: product

const T = Float64

const suite = BenchmarkGroup()
srand(1)

suite["conversions"] = BenchmarkGroup()
rotationtypes = (RotMatrix3{T}, Quat{T}, SPQuat{T}, AngleAxis{T}, RodriguesVec{T})
for (from, to) in product(rotationtypes, rotationtypes)
    if from != to
        name = "$(string(from)) -> $(string(to))"
        # use eval here because of https://github.com/JuliaCI/BenchmarkTools.jl/issues/50#issuecomment-318673288
        suite["conversions"][name] = eval(:(@benchmarkable convert($to, rot) setup = rot = rand($from)))
    end
end

paramspath = joinpath(dirname(@__FILE__), "benchmarkparams.jld")
if isfile(paramspath)
    loadparams!(suite, BenchmarkTools.load(paramspath, "suite"), :evals);
else
    tune!(suite, verbose = true)
    BenchmarkTools.save(paramspath, "suite", params(suite))
end

results = run(suite, verbose=true)
for result in results["conversions"]
    println("$(first(result)):")
    display(minimum(last(result)))
    println()
end
