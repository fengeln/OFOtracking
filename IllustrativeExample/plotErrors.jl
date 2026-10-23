using JLD2
using Plots
using LinearAlgebra
using MAT
pyplot()
Plots.closeall()

colorvec = ["#dc267f", "#fe6100", "#785ef0", "#ffb000", "#648fff"]

datasets1 = []
datasets2 = []
ηvec = 0.01:0.01:1.0  # Step sizes

for η in ηvec
    filename = "data/OFOdata_eta_$(η)_mode_1.jld2"
    push!(datasets1, JLD2.load(filename))
end
for η in ηvec
    filename = "data/OFOdata_eta_$(η)_mode_2.jld2"
    push!(datasets2, JLD2.load(filename))
end

totalerror1 = zeros(length(ηvec),1)
for (i, data) in enumerate(datasets1)
    error1 = data["error"]
    totalerror1[i] = sum(norm(error1))
end
totalerror2 = zeros(length(ηvec),1)
for (i, data) in enumerate(datasets2)
    error2 = data["error"]
    totalerror2[i] = sum(norm(error2))
end
matwrite("data/error1.mat", Dict("totalerror1" => totalerror1))
matwrite("data/error2.mat", Dict("totalerror2" => totalerror2))
tcontrol = datasets1[1]["tcontrol"]

plt0 = Plots.plot(xlabel = "stepsize η", ylabel = "Σ(error)", size=(600, 200), framestyle = :box)
plot!(plt0, ηvec, totalerror1, label="Mode 1", linewidth = 2, color = colorvec[1], xlims = (0.01, 1),ylims = (0.0, 15))
plot!(plt0, ηvec, totalerror2, label="Mode 2", linewidth = 2, color = colorvec[3], xlims = (0.01, 1),ylims = (0.0, 15))