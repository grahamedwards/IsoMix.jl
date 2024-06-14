

"""

    loglikelihood(x::Number, D <: Datum)

Calculate the relative¹ log-likelihood that `x` is drawn from a distribution `D`, which may be normal ([`Norm`](@ref)), lognormal ([`logNorm`](@ref)), or uniform ([`Unf`](@ref)). Note that log-space constants are dropped from the calculations for efficiency (`metropolis only compares log-ratios`)

---

    loglikelihood(c<:Component, d<:Data)

Calculate the loglikelihood that the species composition/concentrations in `c` were drawn from the corresponding distributions in `d`.

---

    loglikelihood(s<:System , p<:Prior)

Calculate the loglikelihood that the components in `s` were drawn from the corresponding `Data` in `p`.

---

    loglikelihood(m::Vector, d::Measurements)
    loglikelihood(m<:Model , d<:DataSet)

Calculate the loglikelihood that the model compositions/concentrations in `m` were drawn from the corresponding measured distribution(s) in `d`.

"""
loglikelihood(x::Number,D::Norm) = -(x-D.m)*(x-D.m)/(2*D.s*D.s)
loglikelihood(x::Number,D::logNorm) = (lnx = log(x); -lnx-(lnx-D.lm)*(lnx-D.lm) / (2*D.ls*D.ls))
loglikelihood(x::Number,D::Unf) = ifelse(D.a <= x <= D.b, -log(D.b-D.a), -Inf)
loglikelihood(x::Number,D::Constant) = ifelse(x===D.x, 0, -Inf)

function loglikelihood(c::C, d::D) where {C<:Component, D <: Data}
    @assert fieldnames(C) == fieldnames(D)
    ll = 0.0
    @inbounds @simd ivdep for i = fieldnames(C) 
        ll += loglikelihood(getfield(c,i),getfield(d,i))
    end
    ll
end

function loglikelihood(s::S, p::P) where {S<:System, P<:Prior}
    @assert fieldnames(S) == fieldnames(P)
    ll = 0.0
    @inbounds for i = fieldnames(S) # probably batch this at some point.
        ll += loglikelihood(getfield(s,i),getfield(p,i))
    end
    ll
end

function loglikelihood(m::Vector{Float64}, d::Measurements)
    ll = 0.0
    @inbounds for i = eachindex(d)
        @inbounds @simd ivdep for j = eachindex(m)
            lli = loglikelihood(m[j],Norm(d.m[i], d.s[i]))
            ll += ifelse(isnan(lli), 0, lli)
        end
    end
    ll
end

function loglikelihood(m::M, d::D) where {M<:Model, D <:DataSet}
    @assert fieldnames(M) == fieldnames(D)
    ll = 0.0

    @inbounds for i = fieldnames(M) # probably batch this at some point.
        ll += loglikelihood(getfield(m,i),getfield(d,i))
    end
    ll
end



