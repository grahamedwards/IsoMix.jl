"""

    initialguess(p<:Prior)

Returns a [`System`](@ref) instance corresponding to provided [`Prior`](@ref) instance `p`, assigning initial guesses of each component composition using the central tendancy of each [`Datum`](@ref): mean of `Norm`, log-mean of `logNorm`, midpoint of `Unf`, value of `Constant`, and an arbitrary value of 0.5 for `Unconstrained`.

"""
function initialguess(p::P) where {D<:Data, P<:Prior{D}}
    prior = ()
    @inbounds for i = fieldnames(P)
        data = ()
        @inbounds for j = fieldnames(D)
            data = (data..., initialguess(getfield(getfield(p,i),j)))
        end
        prior=(prior..., Component(data...))
    end
    return System(prior...)
end
initialguess(x::Norm) = x.m 
initialguess(x::logNorm) = exp(x.lm)
initialguess(x::Unf) = (x.b + x.a)/2
initialguess(x::Constant) = x.x
initialguess(::Unconstrained) = 0.5



"""

    initialjump(p<:Prior)

Returns a [`System`](@ref) instance corresponding to an initial jumping distribution σ for each component's compositions, given provided [`Prior`](@ref) instance `p`. Jump σ are assigned for each [`Datum`](@ref) as follows, the standar deviation (σ) of `Norm`, the log-σ of `logNorm`, ¼ the range of `Unf`, a value of 0 for `Constant`, and an arbitrary value of 0.1 for `Unconstrained`.

"""
function initialjump(p::P) where {D<:Data, P<:Prior{D}}
    prior = ()
    @inbounds for i = fieldnames(P)
        data = ()
        @inbounds for j = fieldnames(D)
            data = (data..., initialjump(getfield(getfield(p,i),j)))
        end
        prior=(prior..., Component(data...))
    end
    return System(prior...)
end
initialjump(x::Norm) = x.s
initialjump(x::logNorm) = x.ls 
initialjump(x::Unf) = (x.b - x.a)/4
initialjump(x::Constant) = 0
initialjump(::Unconstrained) = 0.1



"""

    jump(s<:System, j<:System; rng)

Given a guess of `System` values `s` and corresponding jumping distribution scales in `j`, randomly perturbs one component of `s` and returns the new guess as well as a tuple containing the jumped [`Component`](@ref)/endmember (`:A`, `:B`, `:C`), the jumped component composition (e.g. `:x`, `:y`, `:cx`), and the jump value.

# Example

    julia> jumpedsystem, t = jump(System(Component(1,2),Component(3,4)), System(Component(.1,.5),Component(.2,.4)), rng=IsoMix.Random.Xoshiro(1)); 

    julia> jumpedsystem
    System2{Component1}(Component1(1.0, 2.349413341845734), Component1(3.0, 4.0))

    julia> t
    (:A, :cx, 0.34941334184573425)

"""
function jump(p::S, j::S ; rng=Random.Xoshiro()) where {C<:Component, S<:System{C}}
    si, ci = rand(rng,fieldnames(S)), rand(rng,fieldnames(C))
    j = getfield(getfield(j,si),ci)*rand(rng)
    update(p, si, ci, getfield(getfield(p,si),ci) + j), (si, ci, j)
end


"""

    update(s<:System, si::Symbol, ci::Symbol, v)

Update the `Component` field `ci` of `System` field `si` of `s` with the value `v`.

# Example

    julia> update(System(Component(1,2),Component(3,4)), :B, :x, 12.)
    System2{Component1}(Component1(1.0, 2.0), Component1(12.0, 4.0))

"""
update
function update(s::System3{Component3}, si::Symbol, ci::Symbol, v::Float64)
    X = s.A
    X = ifelse(si==:B,s.B,X)
    X = ifelse(si==:C,s.C,X)

    X = Component3(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx), ifelse(ci==:y,v,X.y), ifelse(ci==:cy,v,X.cy), ifelse(ci==:z,v,X.z), ifelse(ci==:cz,v,X.cz))
    reassigncomponents(si,X,s.A,s.B,s.C)
end

function update(s::System3{Component2}, si::Symbol, ci::Symbol, v::Float64)
    X = s.A
    X = ifelse(si==:B,s.B,X)
    X = ifelse(si==:C,s.C,X)

    X = Component2(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx), ifelse(ci==:y,v,X.y), ifelse(ci==:cy,v,X.cy))
    reassigncomponents(si,X,s.A,s.B,s.C)
end

function update(s::System3{Component1}, si::Symbol, ci::Symbol, v::Float64)
    X = s.A
    X = ifelse(si==:B,s.B,X)
    X = ifelse(si==:C,s.C,X)

    X = Component1(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx))
    reassigncomponents(si,X,s.A,s.B,s.C)
end

function update(s::System2{Component3}, si::Symbol, ci::Symbol, v::Float64)
    X = ifelse(si==:A,s.A,s.B)
    X = Component3(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx), ifelse(ci==:y,v,X.y), ifelse(ci==:cy,v,X.cy), ifelse(ci==:z,v,X.z), ifelse(ci==:cz,v,X.cz))
    reassigncomponents(si,X,s.A,s.B)
end

function update(s::System2{Component2}, si::Symbol, ci::Symbol, v::Float64)
    X = ifelse(si==:A,s.A,s.B)
    X = Component2(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx), ifelse(ci==:y,v,X.y), ifelse(ci==:cy,v,X.cy))
    reassigncomponents(si,X,s.A,s.B)
end

function update(s::System2{Component1}, si::Symbol, ci::Symbol, v::Float64)
    X = ifelse(si==:A,s.A,s.B)
    X = Component1(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx))
    reassigncomponents(si,X,s.A,s.B)
end


function reassigncomponents(si::Symbol,X::T, A::T, B::T) where T<:Component
    if si==:A System2(X,B)
    elseif si==:B System2(A,X)
    else error("System component/endmember must be A or B")
    end
end
function reassigncomponents(si::Symbol,X::T, A::T, B::T, C::T) where T<:Component
    if si==:A System3(X,B,C)
    elseif si==:B System3(A,X,C)
    elseif si==:C System3(A,B,X)
    else error("System component/endmember must be A, B, or C")
    end
end

