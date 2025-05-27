# Simulating univariate time series

Let's simulate univariate time series data from a power spectrum model using the method of [1995A&A...300..707T](@cite).
First we need to load the `Tonari` package and the `Random` package for generating random numbers.

```@example simulation_univariate
using Tonari
using Random
using Plots
using StatsBase
```

Let's define the power spectrum model, for instance we will use [`SingleBendingPowerLaw`](@ref).

```@example simulation_univariate
psd_model = SingleBendingPowerLaw(1.0, 1.13, 1.2e-1, 4.22)
f = 10.0 .^ range(-3, stop=2, length=1000)
psd = psd_model(f)
plot(f, psd, xscale=:log10, yscale=:log10, xlabel="Frequency", ylabel="Power Spectrum Density", label="Power Spectrum", framestyle=:box)
```
## Regularly sampled data

Now we can simulate the time series data. We first define a [`Simulation`](@ref) struct with the power spectrum model,the total time `T`, the time step `Δt`, the mean count rate `μ`, and the variance `σ²`. We can also add the `S_low` and `S_high` parameters to the simulation object. The `S_low` and `S_high` parameters are factors that extend the grid of frequencies. This is useful to avoid the cyclic effect of the discrete Fourier transform, the resulting time series is extracted from a longer time interval.


```@example simulation_univariate

T, Δt = 504.2, 0.132
simu = Simulation(psd_model, T, Δt, 10, 10)
```

Now we can simulate the time series data by sampling the generative model.

```@example simulation_univariate
rng = MersenneTwister(42)
N = 10
t, x, σ = sample(rng, simu, N, error_size = 0.25)
```

!!! note
    The `split_long` parameter is used to split the long time series into smaller segments to speed up the computation, by default this is enabled. This is useful when we want to simulate a large number of time series.

Finally, we can plot the simulated time series data.

```@example simulation_univariate
scatter(t, x[:, 1], label = nothing, yerr = σ[:, 1], xlabel = "Time (days)", ylabel = "Value", framestyle = :box, ms = 2)
```

We can check that the simulated time series data has the same power spectrum as the model.

```@example simulation_univariate
f,I = periodogram(t,x,apply_end_matching=false,subtract_mean=true)
noise_level = 2Δt*mean(σ.^2)
fx = 10.0 .^ range(-3, stop=0.5, length=1000)
plot(f,I,yscale=:log10,xscale=:log10,xlabel="Frequency (Hz)",ylabel="Power",label="Periodogram",framestyle=:box)
plot!(fx,psd_model(fx),label="Model",linewidth=2)
hline!([noise_level],label="Noise level",linewidth=2,linestyle=:dash)
```

## Irregularly sampled data

We can also simulate irregularly sampled data. We first need to define the time stamps `t`. Let's draw the times from a uniform distribution.

```@example simulation_univariate
t = sort(sample(rng,collect(0:Δt:T),150,replace=false))
```

We now define the [`Simulation`](@ref) object  with the new time stamps.

```@example simulation_univariate
simu = Simulation(psd_model, t, 10, 10)
```

 We can now sample the generative model.
```@example simulation_univariate
t_obs, x, σ = sample(rng, simu, error_size = 0.25)
scatter(t_obs, x, label = nothing, yerr = σ, xlabel = "Time (days)", ylabel = "Value", framestyle = :box, ms = 2)
```


## With custom errorbars

It is possible to give the errorbars to randomise the data. The errorbars are given as a vector of the same length as the time stamps.

```@example simulation_univariate
simu = Simulation(psd_model, t)
σs = 0.1 .+ 0.1 * rand(rng, length(t_obs))

t_obs, x, σ = sample(rng, simu, σₓ=σs)
@assert σ==σs
scatter(t_obs, x, label = nothing, yerr = σs, xlabel = "Time (days)", ylabel = "Value", framestyle = :box, ms = 2)
```
