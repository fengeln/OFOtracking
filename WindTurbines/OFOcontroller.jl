using OrdinaryDiffEq, LinearAlgebra, Zygote

function hsys(u, θ)
    P1 = 291.02571699437334 -2.12381624171503*θ + 0.00393299608904857*θ^2
    P2 = 3285.094179230328 -24.406334334196227*θ + 0.045196920505850335*θ^2
    P3 = -0.0246175643292057 + 9.004534046925536e-5*θ + 4.182360819709135e-9*θ^2
    P4 = -4192.573015338467 + 31.143661456208402*θ - 0.0576734578361727*θ^2
    P5 = 1.783366054096041 - 0.006608207496876236*θ + 1.164383670917924e-8*θ^2
    P6 = -0.40403759046076004 + 0.002998667077174924*θ - 5.5530880678997525e-6*θ^2
    p1 = 1.0565908409898734
    p2 = 22.860379555743638
    p4 = -33.80867823291091
    p6 = -0.001546219869637349
    y = zeros(2)
    y[1] = p1 + p2*u[1] + p4*u[1]^2 + p6*u[2]^2
    y[2] = P1 + P2*u[1] + P3*u[2] + P4*u[1]^2 + P5*u[1]*u[2] + P6*u[2]^2
    return y
end
#
#
#
#
function ∇h(u, θ)
    P1 = 291.02571699437334 -2.12381624171503*θ + 0.00393299608904857*θ^2
    P2 = 3285.094179230328 -24.406334334196227*θ + 0.045196920505850335*θ^2
    P3 = -0.0246175643292057 + 9.004534046925536e-5*θ + 4.182360819709135e-9*θ^2
    P4 = -4192.573015338467 + 31.143661456208402*θ - 0.0576734578361727*θ^2
    P5 = 1.783366054096041 - 0.006608207496876236*θ + 1.164383670917924e-8*θ^2
    P6 = -0.40403759046076004 + 0.002998667077174924*θ - 5.5530880678997525e-6*θ^2
    p1 = 1.0565908409898734
    p2 = 22.860379555743638
    p4 = -33.80867823291091
    p6 = -0.001546219869637349
    dy1du1 = p2 + 2*p4*u[1]
    dy1du2 = 2*p6*u[2]
    dy2du1 = P2 + 2*P4*u[1] + P5*u[2]
    dy2du2 = P3 + P5*u[1] + 2*P6*u[2]
    return [dy1du1 dy1du2; dy2du1 dy2du2]  # Jacobian matrix
end
#
#
# Defines the objective function (negative total power of both turbines)
# that we try to minimize
function Φ(y, u)
    z = -sum(y) + (u[1]-0.333)^2 + 0.001*u[2]^2   # we add a small penalty on yaw angles to avoid too large yaw angles
    return z
end
#
#
#
#
function ΦmatT(y, u)
    dy, du = jacobian((yy, uu) -> Φ(yy, uu), y, u)
    return vcat(du', dy')  # column vector stacking ∂Φ/∂u on top of ∂Φ/∂y
end
#
#
#
#
function buildHessian(θ)
    P1 = 291.02571699437334 -2.12381624171503*θ + 0.00393299608904857*θ^2
    P2 = 3285.094179230328 -24.406334334196227*θ + 0.045196920505850335*θ^2
    P3 = -0.0246175643292057 + 9.004534046925536e-5*θ + 4.182360819709135e-9*θ^2
    P4 = -4192.573015338467 + 31.143661456208402*θ - 0.0576734578361727*θ^2
    P5 = 1.783366054096041 - 0.006608207496876236*θ + 1.164383670917924e-8*θ^2
    P6 = -0.40403759046076004 + 0.002998667077174924*θ - 5.5530880678997525e-6*θ^2
    p1 = 1.0565908409898734
    p2 = 22.860379555743638
    p4 = -33.80867823291091
    p6 = -0.001546219869637349
    H11 = -2(p4 + P4) + 2
    H12 = -P5
    H22 = -2(p6 + P6) + 2*0.001
    H = [H11  H12;
         H12  H22]
    return H
end
#
#
#
# Find true minimizer of the objective function for a given wind direction
function minimizer_phi(θ)
    P1 = 291.02571699437334 -2.12381624171503*θ + 0.00393299608904857*θ^2
    P2 = 3285.094179230328 -24.406334334196227*θ + 0.045196920505850335*θ^2
    P3 = -0.0246175643292057 + 9.004534046925536e-5*θ + 4.182360819709135e-9*θ^2
    P4 = -4192.573015338467 + 31.143661456208402*θ - 0.0576734578361727*θ^2
    P5 = 1.783366054096041 - 0.006608207496876236*θ + 1.164383670917924e-8*θ^2
    P6 = -0.40403759046076004 + 0.002998667077174924*θ - 5.5530880678997525e-6*θ^2
    p1 = 1.0565908409898734
    p2 = 22.860379555743638
    p4 = -33.80867823291091
    p6 = -0.001546219869637349
    H = buildHessian(θ)
    g = [-p2-P2-(2/3); -P3]
    # --- Solve H*u + g = 0 ---
    u_star = -H \ g
    return u_star
end
#
# Input: ucurrent ∈ R^1, current control input
# Input: ycurrent ∈ R^2, current measurement of output power
# Output: results ∈ OSQP.Result, result of the quadratic programming problem ???
function calcStep(ucurrent, ycurrent, windcurrent)
    ∇hucurrent = ∇h(ucurrent, windcurrent)
    HmatT = [I transpose(∇hucurrent)]
    Gmat =  [10 0;0 0.001]
    direction = -inv(Gmat) * HmatT * ΦmatT(ycurrent, ucurrent)
    return direction
end
#
# Input: α ∈ R^1, step size
# Input: ucurrent ∈ R^2, current control input
# Input: ycurrent ∈ R^2, current measurement of output power
# Output: unext
function OFOrun(α, ucurrent, ycurrent, windcurrent) 
    σ = calcStep(ucurrent, ycurrent, windcurrent)
    unext = ucurrent + α * σ
    return unext
end