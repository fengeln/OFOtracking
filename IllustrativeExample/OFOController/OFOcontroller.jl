function ∇Φ(u, θ)
    gradient = [u[1] - θ; 2*u[2] - 4*θ]
    return gradient
end
# Input: α ∈ R^1, step size
# Input: ucurrent ∈ R^2, current control input
# Output: step direction σ ∈ R^2
function calcStep(α, u, θ) 
    Dmat =  [1 0;0 1]
    direction = -inv(Dmat) *∇Φ(u, θ)  # Compute the step direction using the inverse of Gmat and the gradient of Φ
    return direction
end
#
#
# Input: α ∈ R^1, step size
# Input: ucurrent ∈ R^2, current control input
# Output: unext
function OFOrun(α, ucurrent, θcurrent) 
    σ = calcStep(α, ucurrent, θcurrent)
    unext = ucurrent + α * σ
    return unext
end