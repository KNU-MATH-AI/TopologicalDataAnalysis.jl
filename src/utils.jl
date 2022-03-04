function noisy_circle_data(N, σ)
    r = 1 .+ σ*randn(1,N)
    θ = rand(Uniform(0, 2π), (1,N))

    rcosθ = r .* cos.(θ)
    rsinθ = r .* sin.(θ)

    data = [rcosθ; rsinθ]

    return data
end