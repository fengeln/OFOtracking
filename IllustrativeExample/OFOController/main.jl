using JLD2, LinearAlgebra
include("OFOcontroller.jl")

function minimizer_phi(theta)
    u_star = [theta;2*theta]
    return u_star
end
function calcRadius(η,λ_min,λ_max,δmax)
    ρ = max(abs(1-η*λ_min), abs(1-η*λ_max))
    r = δmax / (1-ρ)
    return r
end
function getθ(i, mode)
    if mode == 1
        return 0.01*sin(0.9*π*i) + 0.1*cos(0.01*π*i)
    else
        return 12.8e-4*sin(0.99*π*i) + 0.64*cos(0.01*π*i)
    end
end
function calcCost(u,θ)
    return 0.5*(u[1]-θ)^2 + (u[2]-2*θ)^2
end

ηvec = 0.01:0.01:1.0 # Step sizes
titer = 10000       # Number of iterations
mode = 1

for η in ηvec
    println("Running OFO with step size η = $η")
    u_ofo = Vector{Vector{Float64}}()
    u_opt = Vector{Vector{Float64}}()
    cost = Vector{Float64}()
    θ = Vector{Float64}()
    error = Vector{Float64}()
    global ucurrent = minimizer_phi(getθ(0,mode))
    push!(u_ofo, ucurrent)  # Store control inputs
    push!(u_opt, ucurrent)
    push!(error, norm(ucurrent - ucurrent))
    push!(θ, getθ(0,mode))
    push!(cost, calcCost(ucurrent, getθ(0,mode)))
    for i in 1:titer
        global θcurrent = getθ(i, mode)
        global costcurrent = calcCost(ucurrent, θcurrent) 
        unext = OFOrun(η, ucurrent, θcurrent)
        
        push!(u_ofo, [ucurrent[1]; ucurrent[2]])
        push!(u_opt, minimizer_phi(getθ(i,mode)))
        push!(error, norm(ucurrent - minimizer_phi(getθ(i,mode))))
        push!(θ, θcurrent)
        push!(cost, costcurrent)
        global ucurrent = unext
        println("Step $i: Control input = $(round.(ucurrent, digits=2))")
    end


    u1 = getindex.(u_ofo, 1)              # extract first component of ucurrent
    u2 = getindex.(u_ofo, 2)              # extract second component
    uopt1 = getindex.(u_opt, 1)          # extract first component of optimal control
    uopt2 = getindex.(u_opt, 2)          # extract second component of optimal control

    λ_min = 1
    λ_max = 2 
    δmax = maximum(norm([uopt1[t]; uopt2[t]] .- [uopt1[t+1]; uopt2[t+1]]) for t in 1:length(uopt1)-1)
    radius = calcRadius(η, λ_min, λ_max, δmax)
    println("Calculated radius without information: $radius")

    totalerror = sum(error)
    println("Total error: $totalerror")

    tcontrol = 0:1:titer


    # ============================================
    # SAVE EVERYTHING NEEDED FOR Plotting
    # ============================================
    JLD2.save("data/OFOdata_eta_$(η)_mode_$(mode).jld2",
        "eta" , η,
        "u1" , u1,
        "u2" , u2,
        "uopt1" , uopt1,
        "uopt2" , uopt2,
        "error" , error,
        "tcontrol" , tcontrol)
end