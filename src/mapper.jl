using Distributions
using Plots

function noise_circle_data(N)
      r = 1 .+ 0.05randn(N)
      θ = rand(Uniform(0, 2π), N)

      data = zeros(2, N)
      for i ∈ 1:N
            data[:, i] = [r[i]*cos(θ[i]), r[i]*sin(θ[i])]
      end

      return data
end

X = noise_circle_data(200)
scatter(X[1,:], X[2,:])

function filter_1d(X)
      min_xAxes = minimum(X[1, :])
      X[1, :] .- min_xAxes
end

function mapper(X, filter, interval)#, interval, overlap)
      # X : Data
      # filter : filter function
      # interval : interval number of range of filter function
      # overlap : 

      Y = filter(X)

      max_Y = maximum(Y)
      min_Y = minimum(Y)

      for i ∈ 1:interval
            interval

      return max_Y, min_Y
end


mapper(X, filter_1d)