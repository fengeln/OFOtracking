include("OFOcontroller.jl")
using CSV
using DataFrames
using JLD2


df = CSV.read("data/wind_1s.csv", DataFrame)
wind_1s = df.wind
wind_1min = wind_1s[1:60:end]

function calcLambdaMinMax(D)
    λ_min = Inf
    λ_max = -Inf
    for wd in 255:285
        H = buildHessian(wd)
        λ = eigvals(D\H)
        λ_min = min(λ_min, minimum(λ))
        λ_max = max(λ_max, maximum(λ))
    end
    return λ_min, λ_max
end
function calcRadius(η,λ_min,λ_max,δmax)
    ρ = max(abs(1-η*λ_min), abs(1-η*λ_max))
    r = δmax / (1-ρ)
    return r
end

α = 0.1 # Step sizes
D = [10 0;0 0.0004]
titer = 100       # Number of iterations
global tstep = 60

u_ofo = Vector{Vector{Float64}}()
u_opt = Vector{Vector{Float64}}()
windDir = Vector{Float64}()

windcurrent = wind_1min[1]  # initial wind direction
ustart = minimizer_phi(windcurrent)  # initial control input (induction factor, yaw offset)
ucurrent = ustart
push!(u_ofo, [ucurrent[1]; ucurrent[2]])  # Store control inputs
push!(windDir, windcurrent)
ycurrent = hsys(ucurrent, windcurrent)

for i in 1:titer
    global windcurrent = wind_1min[i]
    global ucurrent = OFOrun(α, ucurrent, ycurrent, windcurrent)
    global ycurrent = hsys(ucurrent, windcurrent)
    push!(u_ofo, [ucurrent[1]; ucurrent[2]])
    push!(u_opt, minimizer_phi(windcurrent))
    push!(windDir, windcurrent)
    total_power = sum(ycurrent)
    println("Step $i: Control input = $ucurrent, " * "Power = $ycurrent, " * "Total Power = $total_power, " )
end


u1 = getindex.(u_ofo, 1)              # extract first component of ucurrent
u2 = getindex.(u_ofo, 2)              # extract second component
uopt1 = getindex.(u_opt, 1)          # extract first component of optimal control
uopt2 = getindex.(u_opt, 2)          # extract second component of optimal control

λ_min, λ_max = calcLambdaMinMax(D)
println("λ_min = $λ_min, λ_max = $λ_max")
δmax = maximum(norm([uopt1[t]; uopt2[t]] .- [uopt1[t+1]; uopt2[t+1]]) for t in 1:length(uopt1)-1)
radius = calcRadius(α, λ_min, λ_max, δmax)
println("Calculated radius: $radius")


tsec = 0:1:tstep*(titer)  # time vector in seconds for plotting
# ===================+========================
# SAVE EVERYTHING NEEDED FOR POST-PROCESSING
# ============================================
JLD2.save("data/simulation_OFO_alpha_$(α).jld2",
    "alpha" , α,
    "windDir" , windDir,
    "u1" , u1,
    "u2" , u2,
    "uopt1" , uopt1,
    "uopt2" , uopt2,
    "tsec" , tsec)