function get_theme()
    tw = 1.85
    ts = 10
    mts = 5
    ls = 15
    fontsize_theme = Theme(
        fontsize = 16,
        Axis = (
            yticksmirrored = true,
            xticksmirrored = true,
            xminorticksvisible = true,
            yminorticksvisible = true,
            xminorgridvisible = false,
            xgridvisible = false,
            yminorgridvisible = false,
            ygridvisible = false,
            xminortickalign = 1,
            yminortickalign = 1,
            xtickalign = 1,
            ytickalign = 1, xticklabelsize = ls, yticklabelsize = ls,
            xticksize = ts, xtickwidth = tw,
            xminorticksize = mts, xminortickwidth = tw,
            yticksize = ts, ytickwidth = tw,
            yminorticksize = mts, yminortickwidth = tw,
        ), Lines = (
            linewidth = 2.5,
        )
    )
    return fontsize_theme
end

"""
Plot the light curves in raw counts for the source and background(rescaled by area of course!) and the fractional exposure
"""
function plot_light_curve_rawcounts(dm_filtered, en, prefix = "", suffix = ""; extra_data = nothing, xlimits = nothing)
    instruments = dims(dm_filtered, :instrument)
    fig = Figure(size = (900, 350 * length(instruments)))
    for (i, instr) in enumerate(instruments)
        gi = fig[i, 1] = GridLayout()
        axlc = Axis(gi[1, 1], ylabel = "Net counts", xticklabelsvisible = false, title = String(instr), xminorticks = IntervalsBetween(10), yminorticks = IntervalsBetween(5))
        axdist = Axis(gi[1, 2])
        hidespines!(axdist)
        hidedecorations!(axdist)

        axfracexp = Axis(gi[2, 1], xlabel = "Time (s)", ylabel = "Frac. Exp.", xminorticks = IntervalsBetween(10), xticksmirrored = false)
        ss = CairoMakie.stairs!(axlc, dm_filtered.net[Erange = At(en), instrument = At(instr)], step = :post)
        sb = CairoMakie.stairs!(axlc, dm_filtered.bkg[Erange = At(en), instrument = At(instr)], step = :post)
        fs = CairoMakie.stairs!(axfracexp, dm_filtered.frac[fracexp = At(:srcfracexp), instrument = At(instr)], color = "blue", step = :post)
        fb = CairoMakie.stairs!(axfracexp, dm_filtered.frac[fracexp = At(:bkgfracexp), instrument = At(instr)], color = "green", step = :post)
        d = dm_filtered.net[Erange = At(en), instrument = At(instr)].data
        CairoMakie.hist!(axdist, d[.!ismissing.(d)], direction = :x)

        if !isnothing(xlimits)
            CairoMakie.xlims!(axfracexp, xlimits...)
        end
        if !isnothing(extra_data)
            CairoMakie.vspan!(axlc, extra_data[String(instr)]["src_BTI"][:, 1] .- extra_data["time"][1], extra_data[String(instr)]["src_BTI"][:, 2] .- extra_data["time"][1], alpha = 0.2)
        end
        if i == 1
            leg = Legend(gi[2, 2], [ss, sb, fs, fb], ["net", "bkg", "FRACEXP(src)", "FRACEXP(bkg)"], framevisible = false)
            leg.tellheight = false
            leg.tellwidth = false
        end
        # share x and y for the subplots
        linkyaxes!(axlc, axdist)
        linkxaxes!(axlc, axfracexp)

        # scale histogram
        colsize!(gi, 1, Relative(4 / 5))
        rowsize!(gi, 1, Relative(2 / 3))
        colgap!(gi, Relative(0.0))
        rowgap!(gi, Relative(0.05))

        titlelayout = GridLayout(fig[0, 1], halign = :left, tellwidth = false)
        Label(titlelayout[1, 1], "light curves " * en * " keV" * suffix, halign = :left, fontsize = 22, font = "TeX Bold Makie")
        rowgap!(titlelayout, 0)
    end
    if prefix != ""
        prefix = prefix * "_"
    end
    if suffix != ""
        suffix = "_" * suffix
    end
    save(prefix * "raw_light_curve_" * en * ".pdf", fig)
    save(prefix * "raw_light_curve_" * en * ".png", fig)

    return fig
end


function plot_combined_lightcurves(dm_filtered, prefix = "", suffix = "")

    t = collect(dims(dm_filtered, :time))[:].data
    Δt = t[2] - t[1]
    energies = collect(dims(dm_filtered, :Erange)).data
    netcr = sum(dm_filtered.net, dims = :instrument) / Δt
    err_netcr = sqrt.(sum(dm_filtered.varnet, dims = :instrument)) / Δt
    bkgcr = sum(dm_filtered.bkg, dims = :instrument) / Δt
    err_bkgcr = sqrt.(sum(dm_filtered.varbkg, dims = :instrument)) / Δt
    fig = Figure(size = (1200, 300 * length(energies)))
    axes = []

    for (i, en) in enumerate(dims(dm_filtered, :Erange))

        axis = Axis(fig[i, 1], xlabel = "Time (s)", ylabel = "Rate (Count/s)", title = String(en) * " keV", xminorticks = IntervalsBetween(10))
        push!(axes, axis)

        CairoMakie.errorbars!(axis, t, netcr[Erange = At(en)].data[:], err_netcr[Erange = At(en)].data[:], color = :royalblue4)
        CairoMakie.scatter!(t, netcr[Erange = At(en)].data[:], markersize = 4, strokecolor = :royalblue4, strokewidth = 1.0, color = :white)


        CairoMakie.errorbars!(axis, t, bkgcr[Erange = At(en)].data[:], err_bkgcr[Erange = At(en)].data[:], color = :firebrick3)
        CairoMakie.scatter!(t, bkgcr[Erange = At(en)].data[:], markersize = 4, strokecolor = :firebrick3, strokewidth = 1.0, color = :white)
        if i > 1
            linkxaxes!(axis, axes[i - 1])
        end

    end
    if prefix != ""
        prefix = prefix * "_"
    end
    if suffix != ""
        suffix = "_" * suffix
    end
    # CairoMakie.xlims!(axis,7.515935e8,7.51594e8)
    save(prefix * "combined_light_curves" * suffix * ".pdf", fig)
    save(prefix * "combined_light_curves" * suffix * ".png", fig)
    return fig
end
