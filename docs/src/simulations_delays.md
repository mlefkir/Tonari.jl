# Simulating delayed time series

Let's simulate two time series with a delay

```@example simulation_delay
using Tonari
using Random
using Plots
using StatsBase
rng = MersenneTwister(42)
```

Let's define the individual power spectral densities and the phase delay.
We will use two bending power laws and a constant time delay.

```@example simulation_delay
p1 = SingleBendingPowerLaw(1.0, 0.63, 10e-2, 3.22)
p2 = SingleBendingPowerLaw(1.0, 0.78, 13e-2, 3.81)
Δϕ = ConstantTimeLag(12.35)
```

The cross-spectral density is defined by the two power spectral densities and the phase delay.
```@example simulation_delay
cs = CrossSpectralDensity(p1, p2, Δϕ)
```
We can now simulate the time series using the `Simulation` type.
```@example simulation_delay
T, Δt = 504.2, 0.133
simu = Simulation(cs, T, Δt)
```
and sample 50 time series.
```@example simulation_delay
t, y, yerr = sample(rng,simu,50,error_size=0.15)
```

```@example simulation_delay
scatter(t,y[1][:,1],yerr=yerr[1][:,1],label="1", framestyle=:box)
scatter!(t,y[2][:,1],yerr=yerr[2][:,1],label="2", framestyle=:box, xlabel="Time (days)", ylabel="Value")
```
The cross periodogram can be computed using the `cross_periodogram` function. γ² is the coherence, Δφ the phase difference and Δτ the time delay.
```@example simulation_delay
f, γ², γ²_corrected, Δφ, γ²_err, γ²_corrected_err, Δφ_err, Δτ, Δτ_err, P̄₁, P̄₂, N₁, N₂, n = cross_periodogram(t, y[1], y[2], yerr[1],yerr[2])
```
The individual periodograms can be plotted:
```@example simulation_delay
plot(f,P̄₁,label="P₁")
plot!(f,P̄₂,xscale=:log10,yscale=:log10,label="P₂",xlabel="Frequency (d^-1)",ylabel="Power")
hline!([N₁,N₂],linestyle=:dash,label="noise level", framestyle=:box)
```

The coherence and the phase difference can be plotted as a function of the frequency:
```@example simulation_delay
l = @layout [a; b]
p1 = scatter(f, γ²_corrected, yerr=γ²_corrected_err, label="γ² corrected",marker=:square,ms=4,xscale=:log10,size=(800,400))
p1 = scatter!(f, γ², yerr=γ²_err, label="γ²", xscale=:log10,ylims=(0,1.1),ms=4,ylabel="Coherence",link=:x, framestyle=:box)
p2 = scatter(f, Δτ, yerr=Δτ_err, xscale=:log10, ylabel="Time delay (d)",xlabel="Frequency (d^-1)",link=:x, framestyle=:box,label=nothing)
p2 = hline!([Δϕ.Δτ ],color=:black,label="True delay")
plot(p1, p2, layout=l,size=(800,600))
```


## With a constant phaselag

```@example simulation_delay
p1 = SingleBendingPowerLaw(1.0, 0.63, 10e-2, 3.22)
p2 = SingleBendingPowerLaw(1.0, 0.78, 13e-2, 3.81)
Δϕ = ConstantPhaseLag(56.2,5*1/T)

cs = CrossSpectralDensity(p1, p2, Δϕ)
```

```@example simulation_delay
T, Δt = 504.2, 0.133
simu = Simulation(cs, T, Δt)
```

```@example simulation_delay
t, y, yerr = sample(rng,simu,50,error_size=0.15)
```

```@example simulation_delay
scatter(t,y[1][:,1],yerr=yerr[1][:,1],label="1", framestyle=:box)
scatter!(t,y[2][:,1],yerr=yerr[2][:,1],label="2", framestyle=:box, xlabel="Time (days)", ylabel="Value")
```

```@example simulation_delay
f, γ², γ²_corrected, Δφ, γ²_err, γ²_corrected_err, Δφ_err, Δτ, Δτ_err, P̄₁, P̄₂, N₁, N₂, n = cross_periodogram(t, y[1], y[2], yerr[1],yerr[2])
```

```@example simulation_delay
plot(f,P̄₁,label="P₁")
plot!(f,P̄₂,xscale=:log10,yscale=:log10,label="P₂",xlabel="Frequency (d^-1)",ylabel="Power")
hline!([N₁,N₂],linestyle=:dash,label="noise level", framestyle=:box)
```


```@example simulation_delay
l = @layout [a; b ; c]
p1 = scatter(f, γ²_corrected, yerr=γ²_corrected_err, label="γ² corrected",marker=:square,ms=4,xscale=:log10,size=(800,400))
p1 = scatter!(f, γ², yerr=γ²_err, label="γ²", xscale=:log10,ylims=(0,1.1),ms=4,ylabel="Coherence",link=:x, framestyle=:box)

p2 = scatter(f, Δφ, yerr=Δφ_err, xscale=:log10, ylabel="Time delay (d)",xlabel="Frequency (d^-1)",link=:x, framestyle=:box,label=nothing)
p2 = plot!(f,evaluate(Δϕ,f) .* 2 .*pi .* f )
println(evaluate(Δϕ,f[1]) .* 2 .*pi .* f[1] )

p3 = scatter(f, Δτ, yerr=Δτ_err, xscale=:log10, ylabel="Time delay (d)",xlabel="Frequency (d^-1)",link=:x, framestyle=:box,label=nothing)
p3 = plot!(f,-evaluate(Δϕ,f) )

plot(p1, p2,p3, layout=l,size=(800,900))
```