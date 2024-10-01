# Periodogram computation

The periodogram is an estimator for the power spectral density of a time series. It is defined as the squared magnitude of the discrete Fourier transform of the time series. We use the [FFTW.jl](https://github.com/JuliaMath/FFTW.jl) package to compute the discrete Fourier transform efficiently:

```math
X(f) = \sum_{n=1}^{N} x_n e^{-2\pi i f n \Delta t}
```

where $x_n$ is the time series , $f$ is the frequency, and $\Delta t$ is the time interval between data points. The periodogram is then defined as 

```math
I(f) = \frac{2}{T} |X(f)|^2
```

where $T=N\Delta t$ is the total duration of the time series. The periodogram is an unbiased estimator of the power spectral density.

```@example periodogram
using Tonari
using Random
using Plots
using StatsBase 
rng = MersenneTwister(1234)
```

Let's start by simulating time series data with a single bending power law power spectral density model.
```@example periodogram
psd_model = SingleBendingPowerLaw(1.0, 1.13, 1.2e-1, 4.22)
f = 10.0 .^ range(-3, stop=2, length=1000)
psd = psd_model(f)
T, Δt = 504.2, 0.12
simu = Simulation(psd_model, T, Δt)
t, x, σ = sample(rng, simu, 10, error_size = 0.25)
```

We compute the average periodogram of the time series data using the function [`periodogram`](@ref). It is common to subtract the mean of the time series before computing the periodogram. The function also allows for applying end-matching, which is useful for reducing spectral leakage when the first and last points of the time series are not close to each other in value. 

```@example periodogram
f,I = periodogram(t,x,apply_end_matching=false)
```

We can compare the periodogram to the model power spectral density. The noise level is given by $2\Delta t \times \text{mean}(\sigma^2)$, where $\sigma$ is the error on the time series data. See Appendix A of [2003MNRAS.345.1271V](@cite) for other normalisations and noise levels.

```@example periodogram
noise_level = 2Δt*mean(σ.^2)

fx = 10.0 .^ range(-3, stop=0.5, length=1000)
plot(f,I,yscale=:log10,xscale=:log10,xlabel="Frequency (Hz)",ylabel="Power",label="Periodogram",framestyle=:box)
plot!(fx,psd_model(fx),label="Model",linewidth=2)
hline!([noise_level],label="Noise level",linewidth=2,linestyle=:dash)
```