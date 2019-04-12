using Statistics
using LsqFit

# Translation variables
@inline function pq(ϕ, c)
    imax = length(ϕ)
    p = zeros(imax)
    q = zeros(imax)
    p[1] = ϕ[1] * cos(c)
    q[1] = ϕ[1] * sin(c)
    for i in 2:imax
        p[i] = p[i-1] + ϕ[i-1]*cos(c * (i-1))
        q[i] = q[i-1] + ϕ[i-1]*sin(c * (i-1))
    end
    return p, q
end

# Mean square displacement
@inline function Mn_c(ϕ, c, ncut)
    p, q = pq(ϕ, c)
    N = length(ϕ) - ncut
    Mn = zeros(ncut)
    for n in 1:ncut
        Mn[n] = mean([(p[j+n] - p[j])^2 + (q[j+n] - q[j])^2 for j in 1:N])
    end
    return Mn
end

# Oscillatory term
function Vosc_c(ϕ, c, ncut)
    Eϕ = mean(ϕ)
    return [Eϕ^2 * (1 - cos(n*c))/(1 - cos(c)) for n in 1:ncut]
end

function Dn_c(ϕ, c, ncut)
    return Mn_c(ϕ, c, ncut) - Vosc_c(ϕ, c, ncut)
end

function Dn_c_tilde(ϕ, c, ncut)
    Dn = Dn_c(ϕ, c, ncut)
    return Dn .- min(Dn...)
end

# The asymptotic growth rate
function Kc(ϕ, c, ncut)
    linear_func = (x, p) -> p[1] .+ p[2] .* x
    fit = curve_fit(linear_func, log.(1:ncut), log.(Dn_c_tilde(ϕ, c, ncut)[1:end] .+ 1e-2), [0, 0.5])
    return fit.param[2]
end

@inline function regressionmethod(ϕ, ncut)
    c_range = 0.01:0.01:2π
    return median([Kc(ϕ, c, ncut) for c in c_range])
end

@inline function correlationmethod(ϕ, ncut)
    ξ = 1:ncut;
    c_range = 0.01:0.01:2π
    K_c = [cor(ξ, Dn_c(ϕ, c, ncut)) for c in c_range]
    return median(K_c[.!isnan.(K_c)])
end
