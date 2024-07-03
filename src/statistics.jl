

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

Note that rather than calculating the ℓ of every model composition relative to the corresponding distribution in 'measurements`, the three highest ℓ among all instances in a `model` are summed for each measurement. This prevents `mixtropolis` from converging on small parameter spaces that only just fit around the `measurements`, while also limiting aliasing problems.

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
    ll = 0
    @inbounds Polyester.@batch for i = eachindex(measurements.x)
        mx, mcx, my, mcy, mz, mcz = measurements.x[i], measurements.cx[i], measurements.y[i], measurements.cy[i], measurements.z[i], measurements.cz[i]
        ll1 = ll2 = ll3 = -Inf
        @inbounds @simd for j = eachindex(model.x)
            llj = nanll(model.x[j], mx) + nanll(model.cx[j], mcx) + nanll(model.y[j], my) + nanll(model.cy[j], mcy) + nanll(model.z[j], mz) + nanll(model.cz[j], mcz) 
            ll1, ll2, ll3 = top3(llj, ll1, ll2, ll3)
        end
        ll += log(exp(ll1) + exp(ll2) + exp(ll3))
    end
    ll
end
function loglikelihood(model::Model2, measurements::Measurements2)
    ll = 0
    @inbounds Polyester.@batch for i = eachindex(measurements.x)
        mx, mcx, my, mcy = measurements.x[i], measurements.cx[i], measurements.y[i], measurements.cy[i]
        ll1 = ll2 = ll3 = -Inf
        @inbounds @simd for j = eachindex(model.x)
            llj = nanll(model.x[j], mx) + nanll(model.cx[j], mcx) + nanll(model.y[j], my) + nanll(model.cy[j], mcy) 
            ll1, ll2, ll3 = top3(llj, ll1, ll2, ll3)
        end
        ll += log(exp(ll1) + exp(ll2) + exp(ll3))
    end
    ll
end
function loglikelihood(model::Model1, measurements::Measurements1)
    ll = 0
    @inbounds Polyester.@batch for i = eachindex(measurements.x)
        mx, mcx = measurements.x[i], measurements.cx[i]
        ll1 = ll2 = ll3 = -Inf
        @inbounds @simd for j = eachindex(model.x)
            llj = nanll(model.x[j], mx) + nanll(model.cx[j], mcx)
            ll1, ll2, ll3 = top3(llj, ll1, ll2, ll3)
        end
        ll += log(exp(ll1) + exp(ll2) + exp(ll3))
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
nanll(x::Float64,::Unconstrained) = 1.



"""

    IsoMix.top3(y, x1, x2, x3)

Returns the three highest inputs in decreasing order, as long as `x1`, `x2`, and `x3` are provided in decreasing order. **This is optimized for speed and a very specific use case, so entering x values out of order will result in inaccurate ranking.**

"""
function top3(y::Number, x1::Number, x2::Number, x3::Number)
    x1, x2, x3 = ifelse(zero(y) > y >= x1, (y, x1, x2), (x1, x2, x3))
    x2, x3  = ifelse(zero(y) > x1 > y >= x2, (y, x2), (x2, x3))
    x3 = ifelse(zero(y) > x2 > y > x3, y, x3)
    return x1, x2, x3 
end 
