

"""

    loglikelihood(x::Number, D <: Datum)

Calculate the relative¹ log-likelihood that `x` is drawn from a distribution `D`, which may be normal ([`Norm`](@ref)), lognormal ([`logNorm`](@ref)), or uniform ([`Unf`](@ref)). Note that log-space constants are dropped from the calculations for efficiency (`metropolis only compares log-ratios`)

---

    loglikelihood(c<:EmDraw, d<:Endmember)

Calculate the loglikelihood that the species composition/concentrations in `c` were drawn from the corresponding distributions in `d`.

---

    loglikelihood(s<:SysDraw , p<:Prior)

Calculate the loglikelihood that the components in `s` were drawn from the corresponding `Endmember` in `p`.

---

    loglikelihood(model<:Model , measurements<:Measurements) 

Calculate the likelihood (ℓ) that the (highest likelihood) `model` compositions were drawn from the distributions in `measurements`. Uses `Polyester.@batch`-based multithreading for faster runtimes. 

Note that rather than calculating the ℓ of every model composition relative to the corresponding distribution in 'measurements`, only the highest ℓ among all instances in `model` and a given instance in `measurements` are summed. This prevents `mixtropolis` from converging on small parameter spaces that only just fit around the `measurements`.

"""
loglikelihood(x::Number,D::Norm) = -(x-D.m)*(x-D.m)/(2*D.s*D.s)
loglikelihood(x::Number,D::logNorm) = (lnx = log(x); -lnx-(lnx-D.lm)*(lnx-D.lm) / (2*D.ls*D.ls))
loglikelihood(x::Number,D::Unf) = ifelse(D.a <= x <= D.b, -log(D.b-D.a), -Inf)
loglikelihood(x::Number,D::Constant) = ifelse(float(x)===D.x, 0., -Inf)
loglikelihood(x::Number,::Unconstrained) = 0.

function loglikelihood(c::C, d::D) where {C<:EmDraw, D <: Endmember}
    @assert fieldnames(C) == fieldnames(D)
    ll = 0.0
    @inbounds @simd ivdep for i = fieldnames(C) 
        ll += loglikelihood(getfield(c,i),getfield(d,i))
    end
    ll
end

function loglikelihood(s::S, p::P) where {S<:SysDraw, P<:Prior}
    @assert fieldnames(S) == fieldnames(P)
    ll = 0.0
    @inbounds for i = fieldnames(S)
        ll += loglikelihood(getfield(s,i),getfield(p,i))
    end
    ll
end

function loglikelihood(model::Model3, measurements::Measurements3)
    ll = -Inf
    @inbounds Polyester.@batch for i = eachindex(measurements.x)
        mx, mcx, my, mcy, mz, mcz = measurements.x[i], measurements.cx[i], measurements.y[i], measurements.cy[i], measurements.z[i], measurements.cz[i]
        llj = -Inf
        @inbounds @simd ivdep for j = eachindex(model.x)
            llj = nanll(model.x[j], mx) + nanll(model.cx[j], mcx) + nanll(model.y[j], my) + nanll(model.cy[j], mcy) + nanll(model.z[j], mz) + nanll(model.cz[j], mcz) 
            llj = ifelse(llj > ll, llj, ll)
        end
        ll += llj
    end
    ll
end
function loglikelihood(model::Model2, measurements::Measurements2)
    ll = 0
    @inbounds Polyester.@batch for i = eachindex(measurements.x)
        mx, mcx, my, mcy = measurements.x[i], measurements.cx[i], measurements.y[i], measurements.cy[i]
        lli = -Inf
        @inbounds @simd ivdep for j = eachindex(model.x)
            llj = nanll(model.x[j], mx) + nanll(model.cx[j], mcx) + nanll(model.y[j], my) + nanll(model.cy[j], mcy) 
            lli = ifelse(llj > lli, llj, lli)
        end
        ll += lli
    end
    ll
end
function loglikelihood(model::Model1, measurements::Measurements1)
    ll = -Inf
    @inbounds Polyester.@batch for i = eachindex(measurements.x)
        mx, mcx = measurements.x[i], measurements.cx[i]
        llj = -Inf
        @inbounds @simd ivdep for j = eachindex(model.x)
            llj = nanll(model.x[j], mx) + nanll(model.cx[j], mcx)
            llj = ifelse(llj > ll, llj, ll)
        end
        ll += llj
    end
    ll
end

"""
    IsoMix.nanll(x,D<:Union{Norm,Unconstrained})
Same as [`loglikelihood`](@ref), but returns `0.0` instead of `NaN`.
"""
function nanll(x::Float64,D::Norm)
    y = loglikelihood(x,D)
    return ifelse(isnan(y),0.,y)
end
nanll(x::Float64,::Unconstrained) = 0.
