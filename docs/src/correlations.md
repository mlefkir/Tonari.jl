# Lag and correlations between two time series

The lag and correlation functions are useful for understanding the relationship between two time series. Here we show a method to compute the lag and correlation functions using the interpolated cross-correlation function (ICCF).


## Interpolated cross-correlation function (ICCF)

### Explanation

Add explanation of the ICCF from [1986ApJ...305..175G](@cite)
[1998PASP..110..660P](@cite),
[2004ApJ...613..682P](@cite)
Implementation of `sour` for  [2017ApJ...840...41E](@cite) in R by Simon Vaughan [https://github.com/SimonVaughanDataAndCode/sour/](https://github.com/SimonVaughanDataAndCode/sour/)

### Example
Let's simulate a pair of time series with a known cross-spectral density and compute the ICCF. The ICCF is a method to estimate the cross-correlation function between two time series with different sampling rates. The ICCF is computed by interpolating the cross-correlation function between the two time series.

```@example iccf
using Plots,Tonari,Random

rng = MersenneTwister(4)

p1 = SingleBendingPowerLaw(1.0, 0.63, 10.0e-2, 3.22)
p2 = SingleBendingPowerLaw(1.0, 0.78, 13.0e-2, 4.81)
Δϕ = ConstantTimeLag(18.5)
cs = CrossSpectralDensity(p1, p2, Δϕ)

T, Δt = 404.2, 0.1
simu = Simulation(cs, T, Δt)
t, y, yerr = sample(rng, simu, 1, error_size = 0.05)
N1, N2 = 511, 680
# sample indices
p1 = sort(sample(rng, 1:length(t), N1, replace = false))
p2 = sort(sample(rng, 1:length(t), N2, replace = false))

t₁, y₁, σ₁ = t[p1], y[1][p1, 1], yerr[1][p1, 1]
t₂, y₂, σ₂ = t[p2], y[2][p2, 1], yerr[2][p2, 1]

plot(t₁, y₁, yerr=σ₁, label="Time series 1", xlabel="Time", ylabel="Flux", title="Time series 1 and 2")
plot!(t₂, y₂, yerr=σ₂, label="Time series 2")
```

We can compute the ICCF using the `cross_correlate` function. The function returns the cross-correlation function, the time lags, and the centroid of the cross-correlation function from Monte Carlo simulations of the input time series.

```@example iccf
τ_list, r, q = cross_correlate( t₁, y₁, t₂, y₂, σ₁ = σ₁, σ₂ = σ₂,max_lag=150,Δτ=1, compute_errors = true, n_simulations = 1000)
plot(τ_list, r, label="ICCF", xlabel="Time Lag", ylabel="Correlation", title="Interpolated cross-correlation function")
```
We can compute the centroid of the ICCF and plot it as a vertical line.
```@example iccf
m = r .>=0.8*maximum(r)
τ_peak = sum(τ_list[m].*r[m])/sum(r[m])
```

We can plot the distribution of the centroid of the ICCF and use it for error estimation.
```@example iccf
p = histogram(q, bins=50, label="Time Lag ICCF Centroid", xlabel="Time Lag ", ylabel="Frequency",)
p = vline!([τ_peak], label="Lag centroid")
```

## Discrete Correlation function (DCF)

Still need to be implemented
### Explanation


[1988ApJ...333..646E](@cite)