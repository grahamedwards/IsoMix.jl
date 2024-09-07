"""

    IsoMix.initialguess(p<:Prior)

Returns a [`SysDraw`](@ref) instance corresponding to provided [`Prior`](@ref) instance `p`, assigning initial guesses of each component composition using the central tendancy of each [`Datum`](@ref): mean of `Norm`, log-mean of `logNorm`, midpoint of `Unf`, value of `Constant`, and an arbitrary value of 0.5 for `Unconstrained`.

"""
function initialguess(p::P) where {D<:Endmember, P<:Prior{D}}
    prior = ()
    @inbounds for i = fieldnames(P)
        data = ()
        @inbounds for j = fieldnames(D)
            data = (data..., initialguess(getfield(getfield(p,i),j)))
        end
        prior=(prior..., EmDraw(data...))
    end
    return SysDraw(prior...)
end
initialguess(x::Norm) = x.m 
initialguess(x::logNorm) = exp(x.lm)
initialguess(x::Unf) = (x.b + x.a)/2
initialguess(x::Constant) = x.x
initialguess(::Unconstrained) = 0.5



"""

    IsoMix.initialjump(p<:Prior)

Returns a [`SysDraw`](@ref) instance corresponding to an initial jumping distribution σ for each component's compositions, given provided [`Prior`](@ref) instance `p`. Jump σ are assigned for each [`Datum`](@ref) as follows, the standar deviation (σ) of `Norm`, the log-σ of `logNorm`, ¼ the range of `Unf`, a value of 0 for `Constant`, and an arbitrary value of 0.1 for `Unconstrained`.

"""
function initialjump(p::P) where {D<:Endmember, P<:Prior{D}}
    prior = ()
    @inbounds for i = fieldnames(P)
        data = ()
        @inbounds for j = fieldnames(D)
            data = (data..., initialjump(getfield(getfield(p,i),j)))
        end
        prior=(prior..., EmDraw(data...))
    end
    return SysDraw(prior...)
end
initialjump(x::Norm) = x.s
initialjump(x::logNorm) = x.ls 
initialjump(x::Unf) = (x.b - x.a)/4
initialjump(x::Constant) = 0
initialjump(::Unconstrained) = 0.1



"""

    IsoMix.jump(s<:SysDraw, j<:SysDraw; rng)

Given a guess of `SysDraw` values `s` and corresponding jumping distribution scales in `j`, randomly perturbs one component of `s` and returns the new guess as well as a tuple containing the jumped [`EmDraw`](@ref)/endmember (`:A`, `:B`, `:C`), the jumped component composition (e.g. `:x`, `:y`, `:cx`), and the jump value.

# Example

    julia> jumpedsystem, t = jump(SysDraw(EmDraw(1,2),EmDraw(3,4)), SysDraw(EmDraw(.1,.5),EmDraw(.2,.4)), rng=IsoMix.Random.Xoshiro(1)); 

    julia> jumpedsystem
    SysDraw2{EmDraw1}(EmDraw1(1.0, 2.349413341845734), EmDraw1(3.0, 4.0))

    julia> t
    (:A, :cx, 0.34941334184573425)

"""
function jump(p::S, j::S ; rng=Random.Xoshiro()) where {C<:EmDraw, S<:SysDraw{C}}
    si, ci = rand(rng,fieldnames(S)), rand(rng,fieldnames(C))
    j = getfield(getfield(j,si),ci)*randn(rng)
    update(p, si, ci, getfield(getfield(p,si),ci) + j), (si, ci, abs(j))
end


"""

    update(s<:SysDraw, si::Symbol, ci::Symbol, v)

Update the `EmDraw` field `ci` of `SysDraw` field `si` of `s` with the value `v`.

# Example

    julia> update(SysDraw(EmDraw(1,2),EmDraw(3,4)), :B, :x, 12.)
    SysDraw2{EmDraw1}(EmDraw1(1.0, 2.0), EmDraw1(12.0, 4.0))

"""
update
function update(s::SysDraw3{EmDraw3}, si::Symbol, ci::Symbol, v::Float64)
    X = s.A
    X = ifelse(si==:B,s.B,X)
    X = ifelse(si==:C,s.C,X)

    X = EmDraw3(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx), ifelse(ci==:y,v,X.y), ifelse(ci==:cy,v,X.cy), ifelse(ci==:z,v,X.z), ifelse(ci==:cz,v,X.cz))
    reassigncomponents(si,X,s.A,s.B,s.C)
end

function update(s::SysDraw3{EmDraw2}, si::Symbol, ci::Symbol, v::Float64)
    X = s.A
    X = ifelse(si==:B,s.B,X)
    X = ifelse(si==:C,s.C,X)

    X = EmDraw2(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx), ifelse(ci==:y,v,X.y), ifelse(ci==:cy,v,X.cy))
    reassigncomponents(si,X,s.A,s.B,s.C)
end

function update(s::SysDraw3{EmDraw1}, si::Symbol, ci::Symbol, v::Float64)
    X = s.A
    X = ifelse(si==:B,s.B,X)
    X = ifelse(si==:C,s.C,X)

    X = EmDraw1(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx))
    reassigncomponents(si,X,s.A,s.B,s.C)
end

function update(s::SysDraw2{EmDraw3}, si::Symbol, ci::Symbol, v::Float64)
    X = ifelse(si==:A,s.A,s.B)
    X = EmDraw3(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx), ifelse(ci==:y,v,X.y), ifelse(ci==:cy,v,X.cy), ifelse(ci==:z,v,X.z), ifelse(ci==:cz,v,X.cz))
    reassigncomponents(si,X,s.A,s.B)
end

function update(s::SysDraw2{EmDraw2}, si::Symbol, ci::Symbol, v::Float64)
    X = ifelse(si==:A,s.A,s.B)
    X = EmDraw2(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx), ifelse(ci==:y,v,X.y), ifelse(ci==:cy,v,X.cy))
    reassigncomponents(si,X,s.A,s.B)
end

function update(s::SysDraw2{EmDraw1}, si::Symbol, ci::Symbol, v::Float64)
    X = ifelse(si==:A,s.A,s.B)
    X = EmDraw1(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx))
    reassigncomponents(si,X,s.A,s.B)
end



"""

    IsoMix.reassigncomponents(si::Symbol, X, A, B)
    IsoMix.reassigncomponents(si::Symbol, X, A, B,C)

Returns a `SysDraw` containing `EmDraw`s A and B (and C), replacing the field denoted by `si` with `X`

# Examples
    julia> IsoMix.reassigncomponents(:A, EmDraw(2,2), EmDraw(1,1), EmDraw(1,1))

    SysDraw2{EmDraw1}(EmDraw1(2.0, 2.0), EmDraw1(1.0, 1.0))

"""
function reassigncomponents(si::Symbol,X::T, A::T, B::T) where T<:EmDraw
    if si==:A SysDraw2(X,B)
    elseif si==:B SysDraw2(A,X)
    else error("SysDraw component/endmember must be A or B")
    end
end
function reassigncomponents(si::Symbol,X::T, A::T, B::T, C::T) where T<:EmDraw
    if si==:A SysDraw3(X,B,C)
    elseif si==:B SysDraw3(A,X,C)
    elseif si==:C SysDraw3(A,B,X)
    else error("SysDraw component/endmember must be A, B, or C")
    end
end


"""

    IsoMix.extractsystem(s)

Returns all `EmDraw` field values within the `SysDraw` fields of `s`.

see also: [`IsoMix.extractcomponents`](@ref)

"""
extractsystem(s::SysDraw3) = (extractcomponents(s.A)..., extractcomponents(s.B)..., extractcomponents(s.C)...)
extractsystem(s::SysDraw2) = (extractcomponents(s.A)..., extractcomponents(s.B)...)


"""

    IsoMix.extractcomponents(c)

Return all values from `EmDraw` instance `c`.

see also: [`IsoMix.extractsystem`](@ref)

"""
extractcomponents(c::EmDraw1) = (c.x, c.cx)
extractcomponents(c::EmDraw2) = (c.x, c.cx, c.y, c.cy)
extractcomponents(c::EmDraw3) = (c.x, c.cx, c.y, c.cy, c.z, c.cz)

"""

    IsoMix.extractfields(s)

Return all `SysDraw` field names (e.g. `:A`, `:B`) and all `EmDraw` subfield names (e.g. `:x`, `:cx`) as single `Symbol`s defining each, e.g. `:Ax`, `:Bcz`.

"""
function extractfields(::S) where {C<:EmDraw, S<:SysDraw{C}}
    x = ()
    @inbounds for i in fieldnames(S)
        x= (x..., Symbol.(i,fieldnames(C))...)
    end
    x
end

"""

    IsoMix.stopwatch(i, n, t)

Convenience function for [`porewatermetropolis`] that returns a String reporting the progress at step `i` for total steps `n` with start time `t` (in s since the epoch).

"""
function stopwatch(i::Integer,n::Integer,t::Number)
    pd = 10 * i ÷ n
    bar = "■"^pd * "□"^(10-pd)
    tt = round((time() - t) / 60,digits=2)
    string("0% |", bar,"| 100%  ||  step: $i / $n  ||  time: $tt m")
end