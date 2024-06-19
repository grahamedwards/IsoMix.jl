function mixtropolis(prior::Prior, dataset::DataSet; nodes::Int=1001, burninsteps::Int=10_000, chainsteps::Int=10_000, mixtures::Union{Tuple{Number,Number},Tuple{Number,Number,Number,Number}}=(0,0), updates::Int = 10, rng::Random.AbstractRNG = Random.Xoshiro())

    jumpscale = 2.9

    burnupdate, chainupdate = (burninsteps, chainsteps) .÷ updates
    
    chains = Matrix{Float64}(undef, fieldcount(typeof(prior))*fieldcount(typeof(prior.A)),chainsteps)
    lldist = Vector{Float64}(undef,chainsteps)
    chainacceptance = falses(chainsteps)

    mixtures = mixtures == (0,0) ? ifelse(prior isa Prior3, (0,1,0,1), (0,1) ) : mixtures

    fraction = Fraction(mixtures..., n=nodes)

    p, j = initialguess(prior), initialjump(prior)
    m = Model(fraction,p)

    ϕ = p
    mix!(m,ϕ,fraction)
    ll = loglikelihood(m,dataset) + loglikelihood(ϕ,prior)

    burninacceptance=0
    clock= time()
    println("Burn-in --- ", stopwatch(0,ifelse(iszero(burninsteps),1,burninsteps),clock)); flush(stdout)

    @inbounds for i = Base.OneTo(burninsteps)

        ϕ, jumpinfo = jump(p, j, rng=rng)
        mix!(m,ϕ,fraction)
        llϕ =  loglikelihood(m,dataset) + loglikelihood(ϕ,prior)

        # Decide to accept or reject the proposal
        if log(rand(rng)) < (llϕ-ll) 
            j = update(j,jumpinfo[1], jumpinfo[2], jumpinfo[3]*jumpscale) # update j
            p, ll = ϕ, llϕ  # update proposal and log-likelihood
            burninacceptance=+1        
        end

        if iszero(i % burnupdate) # Update progress
            println("Burn-in --- ", stopwatch(i,burninsteps,clock))
            flush(stdout)
        end
    end


    println("\n\n$burninsteps burn-in steps complete. ℓ = $ll, acceptance rate= $(100burninacceptance÷ifelse(iszero(burninsteps),1,burninsteps)) %.\n\n"); flush(stdout)

    println("Recording chain --- ", stopwatch(0,chainsteps,clock)); flush(stdout)

    @inbounds for i = Base.OneTo(chainsteps)

        ϕ, jumpinfo = jump(p, j, rng=rng)
        mix!(m,ϕ,fraction)
        llϕ =  loglikelihood(m,dataset) + loglikelihood(ϕ,prior)

        # Decide to accept or reject the proposal
        if log(rand(rng)) < (llϕ-ll) 
            j = update(j,jumpinfo[1], jumpinfo[2], jumpinfo[3]*jumpscale) # update j
            p, ll = ϕ, llϕ  # update proposal and log-likelihood
            chainacceptance[i]=1     
        end
    
        chains[:,i] .= extractsystem(p)
        lldist[i] = ll

        if iszero(i % chainupdate) # Update progress
            println("Recording chain --- ", stopwatch(i,chainsteps,clock))
            flush(stdout)
        end
    end

    println("\n\n")
    outnames = (extractfields(p)..., :ll, :accept)
    outtuple = ( (chains[i,:] for i =axes(chains,1))..., lldist, chainacceptance )

    NamedTuple{outnames}(outtuple)
end