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
initialguess(x::logNorm) = x.lm
initialguess(x::Unf) = (x.b + x.a)/2
intialguess(x::Constant) = x.x
initialguess(::Unconstrained) = 0.5





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
intialjump(x::Constant) = 0
initialjump(::Unconstrained) = 0.1


function jump(p::S, j::S ; rng=Random.Xoshiro()) where S<:System
    si, ci = drawjump(p)
    v = getfield(getfield(p,si),ci) + getfield(getfield(j,si),ci)*rand(rng)
    update(p, si, ci, v), (si, ci, v)
end

drawjump(::S, rng=Random.Xoshiro()) where {C<:Component, S<:System{C}} = rand(rng,fieldnames(S)), rand(rng,fieldnames(C))



function update(s::System3{Component3}, si::Symbol, ci::Symbol, v::Float64)
    X = s.A
    X = ifelse(si==:B,s.B,X)
    X = ifelse(si==:C,s.C,X)

    X = Component3(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx), ifelse(ci==:y,v,X.y), ifelse(ci==:cy,v,X.cy), ifelse(ci==:z,v,X.z), ifelse(ci==:cz,v,X.cz))

    if si==:A System3(X,B,C)
    elseif si==:B System3(A,X,C)
    elseif si==:C System3(A,B,X)
    else error("System component/endmember must be A, B, or C")
    end
end

function update(s::System3{Component2}, si::Symbol, ci::Symbol, v::Float64)
    X = s.A
    X = ifelse(si==:B,s.B,X)
    X = ifelse(si==:C,s.C,X)

    X = Component2(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx), ifelse(ci==:y,v,X.y), ifelse(ci==:cy,v,X.cy))

    if si==:A System3(X,B,C)
    elseif si==:B System3(A,X,C)
    elseif si==:C System3(A,B,X)
    else error("System component/endmember must be A, B, or C")
    end
end

function update(s::System3{Component1}, si::Symbol, ci::Symbol, v::Float64)
    X = s.A
    X = ifelse(si==:B,s.B,X)
    X = ifelse(si==:C,s.C,X)

    X = Component1(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx))

    if si==:A System3(X,B,C)
    elseif si==:B System3(A,X,C)
    elseif si==:C System3(A,B,X)
    else error("System component/endmember must be A, B, or C")
    end
end

function update(s::System2{Component3}, si::Symbol, ci::Symbol, v::Float64)
    X = ifelse(si==:A,s.A,x.B)
    X = Component3(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx), ifelse(ci==:y,v,X.y), ifelse(ci==:cy,v,X.cy), ifelse(ci==:z,v,X.z), ifelse(ci==:cz,v,X.cz))

    if si==:A System2(X,B)
    elseif si==:B System2(A,X)
    else error("System component/endmember must be A or B")
    end
end

function update(s::System2{Component2}, si::Symbol, ci::Symbol, v::Float64)
    X = ifelse(si==:A,s.A,x.B)
    X = Component2(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx), ifelse(ci==:y,v,X.y), ifelse(ci==:cy,v,X.cy))

    if si==:A System2(X,B)
    elseif si==:B System2(A,X)
    else error("System component/endmember must be A or B")
    end
end

function update(s::System2{Component1}, si::Symbol, ci::Symbol, v::Float64)
    X = ifelse(si==:A,s.A,x.B)
    X = Component1(ifelse(ci==:x,v,X.x), ifelse(ci==:cx,v,X.cx))

    if si==:A System2(X,B)
    elseif si==:B System2(A,X)
    else error("System component/endmember must be A or B")
    end
end





