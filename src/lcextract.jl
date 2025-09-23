SUPPORTED_OBSERVATORIES = [:nustar, :xmmnewton]

"""
    get_bad_time_intervals(header_gti, gti_df, t0_curr,tend_curr)

Extract the bad time intervals (BTIs) from the Good Time Intervals (GTIs).

# Arguments
- `header_gti` : FITS header of the Good Time Intervals
- `gti_df` : DataFrame of the Good Time Intervals
- `t0_curr::Float64` : Start time of the time series
- `tend_curr::Float64` : End time of the time series

# Returns
- `BTI::Array{Float64,2}` : Matrix of dimension N*2 where N is the number of intervals. The other dimension contains START and STOP times of the BTIs.
"""
function get_bad_time_intervals(header_gti, gti_df, t0_curr, tend_curr)

    n_gti = header_gti["NAXIS2"]

    # if we have only one GTI
    if n_gti == 1
        if gti_df[1, "START"] < t0_curr # if the GTI starts before the exposure
            if gti_df[1, "STOP"] > tend_curr # if the GTI ends after the end the exposure
                BTI = nothing
            else
                BTI = [[gti_df[1, "STOP"], tend_curr]] # .- t0_curr
            end
        else  # if the gti is between the start and STOP
            ## CHECK THIS LINE AS WE MAY BE MISSING SOMETHING?
            BTI = [[0, gti_df[1, "START"]] [gti_df[1, "STOP"], tend_curr]] # .- t0_curr
        end
    else
        BTI = zeros(n_gti, 2)
        for j in 1:(n_gti - 1)
            BTI[j, 1] = gti_df[j, "STOP"] #- t0_curr
            BTI[j, 2] = gti_df[j + 1, "START"] #- t0_curr
        end
        # add a bad time after the last point to ensure that we don't collect events after the last bin
        BTI[end, 1] = gti_df[end, "STOP"]
        BTI[end, 2] = tend_curr
    end
    return BTI
end

"""
    get_fractional_exposure(t_bin, BTI, Δt)

Compute the fractional exposure for each bin. The fractional exposure is defined as:
    - 1 if the bin is totally good (not in a BTI)
    - 0 if the bin is totally bad (fully in a BTI)
    - r if the bin is partially in a BTI

    where r = good duration of the bin / Δt

# Arguments
- `t_bin::Vector{Float64}` : Time array.
- `BTI::Array{Float64,2}` : BTIs used to compute the
- `Δt::Float64` : Time bin duration

# Returns
- `BTI::Array{Float64,2}` : Matrix of dimension N*2 where N is the number of intervals. The other dimension contains START and STOP times of the BTIs.
"""
function get_fractional_exposure(t_bin, BTI, Δt)
    Frac_EXP = ones(size(t_bin))
    n_bti = size(BTI, 1)
    for i in 1:n_bti
        lower = BTI[i, 1]
        higher = BTI[i, 2]

        bti_start = searchsortedfirst(t_bin, lower) - 1 #
        bti_end = searchsortedfirst(t_bin, higher) - 1 # this is the beginning of the last bti bin! This means that t_bin[bti_end] is also bad!
        bti_start = max(0, bti_start)

        if bti_end >= 0
            badbin_size = bti_end - bti_start
            if badbin_size == 0
                Frac_EXP[bti_start] = Frac_EXP[bti_start] - (higher - lower) / Δt
            end
            if badbin_size > 1
                Frac_EXP[(bti_start + 1):(bti_end - 1)] .= 0.0
            end
            if badbin_size > 0
                Frac_EXP[bti_start] = Frac_EXP[bti_start] - (t_bin[bti_start] + Δt - lower) / Δt
                Frac_EXP[bti_end] = Frac_EXP[bti_end] - (higher - t_bin[bti_end]) / Δt
            end
        end
    end
    return Frac_EXP
end

"""
    get_events(header_gti, hdu_src, hdu_bkg, Δt_user::Float64, PATTERN::Int64, observatory::Symbol; FLAG::Int64=0, t_clip_start::Float64=30., t_clip_stop::Float64=30.,filename_src="",filename_bkg="")

# Arguments
- `header_gti` : FITS header of the Good Time Intervals
- `hdu_src` : FITS hdu for the source event list
- `hdu_bkg` : FITS hdu for the background event list
- `Δt_user::Float64` : input time bin duration in seconds
- `PATTERN::Int64` : maximum pattern
- `observatory:Symbol` : observatory: :xmm or :nustar
- `FLAG:::Int64` : flag to keep for the event list
- `t_clip_start::Float64` : duration to clip at the beginning of the event list
- `t_clip_stop::Float64` : duration to clip at the end of the event list

# Returns
- `bin_list::Vector{Float64}` : Time bin array.
- `t0::Float64` : first event recorded in the event list
- `t_start::Float64` : start time of the bin_list, this accounts for the clipping
- `t_stop::Float64` : stop time of the bin_list, this accounts for the clipping
- `time_src::Vector{Float64}` : List of event arrival times for the source
- `pi_src::Vector{Float64}` : List of event energies for the source
- `time_bkg::Vector{Float64}` : List of event arrival times for the background
- `pi_bkg::Vector{Float64}` : List of event energies for the background
- `Δt::Float64` : time bin duration obtained
- `filename_src::String` : Filename of the source event list
- `filename_bkg::String` : Filename of the background event list
"""
function get_events(hdu_src, hdu_bkg, Δt_user::Float64, PATTERN::Int64, observatory::Symbol, CCDNR_src::Int64, CCDNR_bkg::Int64; FLAG = 0, t_clip_start::Float64 = 0.0, t_clip_stop::Float64 = 0.0, filename_src = "", filename_bkg = "", instr::String = "")


    if observatory == :xmmnewton
        gti = hdu_src[GTI_extension(observatory, CCDNR_format(CCDNR_src))]
        header_gti = read_header(gti)

        FRAME_TIME = header_gti["FRMTIME"] #frame time in miliseconds
        quality_name = "FLAG"
        pattern_name = "PATTERN"
    elseif observatory == :nustar
        FRAME_TIME = 0.1 # miliseconds
        quality_name = "STATUS"
        pattern_name = "GRADE"
    end

    if observatory == :nustar
        @assert filename_src != "" && filename_bkg != "" "The filenames for the src and bkg event lists must be provided if we want to extract the events"
        time_src, pi_src, flag_src, pattern_src, time_bkg, pi_bkg, flag_bkg, pattern_bkg = get_nustar_events(filename_src, filename_bkg)
        println("Total events in source region: $(length(time_src))")
        println("Total events in background region: $(length(time_bkg))")

    else
        time_src, pi_src, flag_src, pattern_src, ccd_src = read.(Ref(hdu_src["EVENTS"]), ["TIME", "PI", quality_name, pattern_name, "CCDNR"])
        time_bkg, pi_bkg, flag_bkg, pattern_bkg, ccd_bkg = read.(Ref(hdu_bkg["EVENTS"]), ["TIME", "PI", quality_name, pattern_name, "CCDNR"])

        # we need to check that the events are on the same CCD
        @assert length(unique(ccd_src)) == 1 && Int64(ccd_src[1]) == CCDNR_src "The source events of $(instr) are not all on the same CCD, expected $(CCDNR_src) but got $(Int.(unique(ccd_src)))"

        @assert length(unique(ccd_bkg)) == 1 && Int64(ccd_bkg[1]) == CCDNR_bkg "The background events of $(instr) are not all on the same CCD, expected $(CCDNR_bkg) but got $(Int.(unique(ccd_bkg)))"
    end


    N_frames_per_bin = floor(Int, Δt_user * 1.0e3 / FRAME_TIME)
    Δt = N_frames_per_bin * FRAME_TIME / 1.0e3 # true time bin duration

    q_src = (pattern_src .<= PATTERN)
    q_bkg = (pattern_bkg .<= PATTERN)
    if observatory != :nustar
        q_src .&= (flag_src .<= FLAG)
        q_bkg .&= (flag_bkg .<= FLAG)
    end

    t0 = time_src[q_src][1]

    time_src, pi_src = time_src[q_src], pi_src[q_src]
    time_bkg, pi_bkg = time_bkg[q_bkg], pi_bkg[q_bkg]

    println("Events in source region after pattern and flag filtering: $(length(time_src))")
    println("Events in background region after pattern and flag filtering: $(length(time_bkg))")

    # t_max = time_src[end]
    t_start = read_header(hdu_src["EVENTS"])["TSTART"] + t_clip_start
    t_stop = read_header(hdu_src["EVENTS"])["TSTOP"] - t_clip_stop + Δt
    bin_list = range(t_start, t_stop, step = Δt)
    return bin_list, t0, t_start, t_stop, time_src, pi_src, time_bkg, pi_bkg, Δt
end

"""
    get_nustar_events(filename_src, filename_bkg)

This function collects NuSTAR events using CFITSIO as the current version of FITSIO
does not allow for BitArrays columns.

# Arguments
- `filename_src::String` : Path to the source event list
- `filename_bkg::String` : Path to the background event list

# Returns
- `time_src::Vector{Float64}` : Recorded event time in the source region
- `pi_src::Vector{Float64}` : Recorded PI in the source
- `flag_src::Vector{Int16}` : Recorded Flag in the source
- `pattern_src::Vector{Int16}` : Recorded pattern in the source
- `time_bkg::Vector{Float64}` : event time in the background region
- `pi_bkg::Vector{Float64}` : Recorded PI in the background
- `flag_bkg::Vector{Int16}` : Recorded Flag in the background
- `pattern_bkg::Vector{Int16}` : Recorded pattern in the background

"""
function get_nustar_events(filename_src::String, filename_bkg::String)

    hdu_src = fits_open_table(filename_src)
    fits_movnam_hdu(hdu_src, "EVENTS")

    n_fields = parse(Int, fits_read_keyword(hdu_src, "TFIELDS")[1])
    n_events = parse(Int, fits_read_keyword(hdu_src, "NAXIS2")[1])

    # read NuSTAR event list with CFITSIO as FITSIO does not support BitArrays yet
    cols_names = ["TIME", "PI", "STATUS", "GRADE"]

    col_indexes = Dict()
    for i in 1:n_fields
        u = replace(replace(fits_read_keyword(hdu_src, "TTYPE$i")[1], " " => ""), "'" => "")
        if u in cols_names
            col_indexes[u] = i
        end
    end

    time_src, pi_src, flag_src, pattern_src = zeros(n_events), zeros(n_events), zeros(Int16, n_events), zeros(Int16, n_events)

    fits_read_col(hdu_src, col_indexes["PI"], 1, 1, pi_src)
    fits_read_col(hdu_src, col_indexes["TIME"], 1, 1, time_src)
    fits_read_col(hdu_src, col_indexes["STATUS"], 1, 1, flag_src)
    fits_read_col(hdu_src, col_indexes["GRADE"], 1, 1, pattern_src)

    hdu_bkg = fits_open_table(filename_bkg)
    fits_movnam_hdu(hdu_bkg, "EVENTS")
    n_events = parse(Int, fits_read_keyword(hdu_bkg, "NAXIS2")[1])

    time_bkg, pi_bkg, flag_bkg, pattern_bkg = zeros(n_events), zeros(n_events), zeros(Int16, n_events), zeros(Int16, n_events)

    fits_read_col(hdu_bkg, col_indexes["PI"], 1, 1, pi_bkg)
    fits_read_col(hdu_bkg, col_indexes["TIME"], 1, 1, time_bkg)
    fits_read_col(hdu_bkg, col_indexes["STATUS"], 1, 1, flag_bkg)
    fits_read_col(hdu_bkg, col_indexes["GRADE"], 1, 1, pattern_bkg)

    return time_src, pi_src, flag_src, pattern_src, time_bkg, pi_bkg, flag_bkg, pattern_bkg
end
"""
    interpolate_randomise_counts(rng, t, counts, mask::BitVector)

Interpolate the points in the time series and randomise with Poisson noise.

# Arguments
- `rng::Random.AbstractRNG`: Random number generator.
- `t::Vector{Float64}` : time array.
- `counts::Vector{Int64}` : counts array.
- `mask::Vector{Bool}` : mask for the points we want to fill in the time series, i.e. 0 for the bad/missing values and 1 for the good values
"""
function interpolate_randomise_counts(rng::Random.AbstractRNG, t, counts::Vector{Int64}, mask::BitVector)
    x_interp = linear_interpolation(t[mask], counts[mask])(t[.!mask])
    return x_interp = rand.(rng, Poisson.(x_interp))
end


"""
    get_energy_array(energies, observatory)

Returns the array containing the PI energies for the light curve extraction for the given observatory.

# Arguments
- `energies::Vector{Float64}` : energies in keV
- `observatory::Symbol` : observatory either :xmmnewton or :nustar

# Returns
- `PI::Vector{Float64}` : PIs for the event list.
"""
function get_energy_array(energies::Vector{Float64}, observatory::Symbol)
    if observatory == :xmmnewton
        PI = energies .* 1.0e3
    elseif observatory == :nustar
        PI = (energies .- 1.6) / 0.04
    else
        error("The observatory $observatory is not supported")
    end
    return PI
end

"""
    get_pairs(PI)
"""
function get_pairs(PI)
    bands = []
    for j in 2:length(PI)
        en = [PI[j - 1], PI[j]]
        push!(bands, en)
    end
    return bands
end
"""
    get_energy_bands(energies)

Generate a list of energy bands and add a band for the min and max energies.

# Returns
- `Vector{String}` : List of energy bands in keV.

# Example
julia> energies = [0.2,1.5,10.0]
julia> get_energy_bands(energies)
["0.2-1.5", "1.5-10.0", "0.2-10.0"]
"""
function get_energy_bands(energies)
    bands = []
    for j in eachindex(energies)
        if j == length(energies)
            en = "$(minimum(energies))-$(maximum(energies))"
        else
            en = "$(energies[j])-$(energies[j + 1])"
        end
        push!(bands, en)
    end
    return bands
end

function CCDNR_format(CCDNR::Int64)
    if CCDNR < 10
        return "0$CCDNR"
    else
        return "$CCDNR"
    end
end

function GTI_extension(observatory, curr_CCDNR)
    if observatory == :xmmnewton
        return "STDGTI$curr_CCDNR"
    elseif observatory == :nustar
        return "STDGTI"
    else
        error("observatory $observatory not supported")
    end
end
"""
    extract_long_tbin(path::String, source::String, instruments::Vector{String}, Δt_user::Float64, default_instr::String,observatory::Symbol,CCDNR_src::Int64)

Extract the list of bin times from the available instruments for the . The minimum and maximum of the bin times are obtained as the earliest time from and maximum.

The name of the event lists should follow the pattern:

    "path/source_instr_evts_src.fits" and
    "path/source_instr_evts_bkg.fits"

# Arguments
- `path::String` : Path to the directory containing the event lists
- `source::String` : Name of the source
- `instruments::Vector{String}` : List of instruments
- `Δt_user::Float64` : Input time bin duration
- `default_instr::String` : Default instrument for the choice of Δt

# Returns
- `bin_list::Vector{Float64}` : Time bin array.
- `t0::Float64` : first event recorded in the event list
- `min_t_start::Float64` : start time of the bin_list, this accounts for the clipping
- `t_stop::Float64` : stop time of the bin_list, this accounts for the clipping
- `time_src::Vector{Float64}` : List of event arrival times for the source
- `pi_src::Vector{Float64}` : List of event energies for the source
- `time_bkg::Vector{Float64}` : List of event arrival times for the background
- `pi_bkg::Vector{Float64}` : List of event energies for the background
- `Δt::Float64` : time bin duration obtained
- `CCDNR_src::Int64` : CCD number for the source region, usually 4 for XMM-Newton on-axis sources
"""
function extract_long_tbin(path::String, source::String, instruments::Vector{String}, Δt_user::Float64, default_instr::String, observatory::Symbol, CCDNR_src::Union{Vector{Int64}, Int64}, CCDNR_bkg::Union{Vector{Int64}, Int64}; t_clip_start::Float64 = 0.0, t_clip_stop::Float64 = 0.0)
    Δt = 0
    t0_src, tend_src = 0, 0
    t0_src_list, tend_src_list, t0_bkg_list, tend_bkg_list = [], [], [], []
    for (k, instr) in enumerate(instruments)

        filename_src = "$(path)/$(source)_$(instr)_evts_src.fits"
        filename_bkg = "$(path)/$(source)_$(instr)_evts_bkg.fits"
        @assert isfile(filename_src) "src filename $(filename_src) is not found"
        @assert isfile(filename_bkg) "bkg filename $(filename_bkg) is not found"

        hdu_src = FITS(filename_src, "r")
        hdu_bkg = FITS(filename_bkg, "r")
        if CCDNR_src isa Vector
            curr_CCDNR_src, curr_CCDNR_bkg = CCDNR_src[k], CCDNR_bkg[k]
        else
            curr_CCDNR_src, curr_CCDNR_bkg = CCDNR_src, CCDNR_bkg
        end
        # if CCDNR_src isa Vector
        #     _, t0, t_start, t_end, _, _, _, _, Δt_i = get_events( hdu_src, hdu_bkg, Δt_user, PATTERN, observatory,CCDNR_src[k],CCDNR_bkg[k],FLAG=0,filename_src=filename_src,filename_bkg=filename_bkg,instr=instr,t_clip_start=t_clip_start,t_clip_stop=t_clip_stop)
        # else
        #     _, t0, t_start, t_end, _, _, _, _, Δt_i = get_events( hdu_src, hdu_bkg, Δt_user, PATTERN, observatory,CCDNR_src,CCDNR_bkg,FLAG=0,filename_src=filename_src,filename_bkg=filename_bkg,instr=instr,t_clip_start=t_clip_start,t_clip_stop=t_clip_stop)
        # end

        # if instr == default_instr
        #     Δt = Δt_i
        # end
        # min_t_start = min(min_t_start, t_start)
        # min_t0 = min(t0, min_t0)
        # max_t_end = max(t_end, max_t_end)

        if instr == default_instr

            if observatory == :xmmnewton
                gti = hdu_src[GTI_extension(observatory, CCDNR_format(curr_CCDNR_src))]
                header_gti = read_header(gti)

                FRAME_TIME = header_gti["FRMTIME"] #frame time in miliseconds
                quality_name = "FLAG"
                pattern_name = "PATTERN"
            elseif observatory == :nustar
                FRAME_TIME = 0.1 # miliseconds
                quality_name = "STATUS"
                pattern_name = "GRADE"
            end


            N_frames_per_bin = floor(Int, Δt_user * 1.0e3 / FRAME_TIME)
            Δt = N_frames_per_bin * FRAME_TIME / 1.0e3 # true time bin duration
        end

        t0_src = read_header(hdu_src["EVENTS"])["TSTART"]
        tend_src = read_header(hdu_src["EVENTS"])["TSTOP"]
        t0_bkg = read_header(hdu_bkg["EVENTS"])["TSTART"]
        tend_bkg = read_header(hdu_bkg["EVENTS"])["TSTOP"]
        push!(t0_src_list, t0_src)
        push!(t0_bkg_list, t0_bkg)
        push!(tend_bkg_list, tend_bkg)
        push!(tend_src_list, tend_src)

    end
    # if length(unique(t0_src_list))

    bin_list = range(t0_src, tend_src + Δt, step = Δt)
    return bin_list, t0_src, tend_src, Δt
end

"""


# Arguments
- `path::String` : Directory where the event lists and spectra are located
- `source::String` : Name of the source
- `obsid::String` : Observation identifier
- `observatory::Symbol` : Name of the observatory, either


"""
function extract_light_curves(
        path::String,
        source::String, obsid::String,
        observatory::Symbol,
        instruments::Vector{String},
        energies::Vector{Float64},
        Δt_user::Float64,
        CCDNR_src::Union{Vector{Int64}, Int64},
        CCDNR_bkg::Union{Vector{Int64}, Int64};
        default_instr = nothing, min_Frac_EXP = 0.4,
        pattern_max = nothing,
        interpolate = true, t_clip_start::Float64 = 0.0,
        t_clip_stop::Float64 = 0.0, energies_overlap = true, rng = Random.GLOBAL_RNG
    )

    @assert observatory ∈ SUPPORTED_OBSERVATORIES "The observatory provided $(observatory) is not in the list of supported observatories $SUPPORTED_OBSERVATORIES"

    @assert all(energies .> 0) "Check your array of energies!"
    PI = get_energy_array(energies, observatory)

    # if we want to include
    flag = 0
    if isnothing(default_instr)
        default_instr = instruments[1]
    end

    if observatory == :nustar
        pattern_name = "GRADE"
        if isnothing(pattern_max)
            pattern_max = 26
        end
    elseif observatory == :xmmnewton
        pattern_name = "PATTERN"
        if isnothing(pattern_max)
            pattern_max = 4
        end
    end

    if CCDNR_src isa Vector
        @assert length(CCDNR_src) == length(CCDNR_bkg) "The length of CCDNR is not the same for source and background"
    end
    # index for the default_instr instr
    # idx_default = findall(x->x == default_instr, instruments)[1]
    # obtain the binning time
    bin_list, min_t_start, max_t_end, Δt = extract_long_tbin(path, source, instruments, Δt_user, default_instr, observatory, CCDNR_src, CCDNR_bkg, t_clip_start = t_clip_start, t_clip_stop = t_clip_stop)


    data = Dict()
    t_bin = 0

    for (k, instr) in enumerate(instruments)
        @info "Working on $instr"
        data[instr] = Dict()

        # load files
        filename_src = "$(path)/$(source)_$(instr)_evts_src.fits"
        filename_bkg = "$(path)/$(source)_$(instr)_evts_bkg.fits"
        isfile(filename_src)
        isfile(filename_bkg)
        if observatory == :xmmnewton
            filename_backscale_src = "$(path)/$(obsid)_$(source)_$(instr)_spectrum_src_BACKSCALE.fits"
            filename_backscale_bkg = "$(path)/$(obsid)_$(source)_$(instr)_spectrum_bkg_BACKSCALE.fits"
        elseif observatory == :nustar
            filename_backscale_src = "$(path)/$(source)_$(instr)_spec_src.fits"
            filename_backscale_bkg = "$(path)/$(source)_$(instr)_spec_bkg.fits"
        end
        isfile(filename_backscale_src)
        isfile(filename_backscale_bkg)

        backscale_src = read_key(FITS(filename_backscale_src)["SPECTRUM"], "BACKSCAL")[1]
        backscale_bkg = read_key(FITS(filename_backscale_bkg)["SPECTRUM"], "BACKSCAL")[1]
        region_scale = backscale_src / backscale_bkg


        hdu_src = FITS(filename_src, "r")
        hdu_bkg = FITS(filename_bkg, "r")

        if CCDNR_src isa Vector
            curr_CCDNR_src, curr_CCDNR_bkg = CCDNR_src[k], CCDNR_bkg[k]
        else
            curr_CCDNR_src, curr_CCDNR_bkg = CCDNR_src, CCDNR_bkg
        end
        _, t0, t_start, t_end, time_src, pi_src, time_bkg, pi_bkg = get_events(hdu_src, hdu_bkg, Δt_user, pattern_max, observatory, curr_CCDNR_src, curr_CCDNR_bkg, FLAG = flag, filename_src = filename_src, filename_bkg = filename_bkg, instr = instr, t_clip_start = t_clip_start, t_clip_stop = t_clip_stop)

        index_observations = t_start .< bin_list .< t_end

        data[instr]["t0"] = t0
        data[instr]["t_start"] = t_start
        data[instr]["t_end"] = t_end

        # t_start, t_end = min_t_start, max_t_end
        if length(energies) == 2 || energies_overlap
            @info "Assuming we can overlap energy bands in light curves"
            h2_src, h2_bkg = [], []
            if length(energies) == 2
                pairs_PI = [PI]
            else
                pairs_PI = get_pairs(PI)
            end
            for j in 1:length(pairs_PI)
                q_src = pairs_PI[j][1] .<= pi_src .<= pairs_PI[j][2]
                q_bkg = pairs_PI[j][1] .<= pi_bkg .<= pairs_PI[j][2]
                time_src_ = time_src[q_src]
                time_bkg_ = time_bkg[q_bkg]
                h1_src = Hist1D(time_src_; binedges = bin_list, counttype = Int, overflow = false)
                h1_bkg = Hist1D(time_bkg_; binedges = bin_list, counttype = Int, overflow = false)
                push!(h2_src, h1_src)
                push!(h2_bkg, h1_bkg)
            end
        else
            @info "Light curves are extracted in non-overlapping energy bands"

            h2_src = Hist2D((time_src, pi_src); binedges = (bin_list, PI), counttype = Int, overflow = false)
            h2_bkg = Hist2D((time_bkg, pi_bkg); binedges = (bin_list, PI), counttype = Int, overflow = false)
        end
        q_src = minimum(PI) .<= pi_src .<= maximum(PI)
        q_bkg = minimum(PI) .<= pi_bkg .<= maximum(PI)

        h_src_full = Hist1D(time_src[q_src]; binedges = bin_list, counttype = Int, overflow = false)
        h_bkg_full = Hist1D(time_bkg[q_bkg]; binedges = bin_list, counttype = Int, overflow = false)

        # for each energy band
        for j in 1:length(PI)
            if j == length(PI)
                counts_src = bincounts(h_src_full)[:]
                counts_bkg = bincounts(h_bkg_full)[:]
                energy = "$(minimum(energies))-$(maximum(energies))"
            else
                if h2_src isa Hist2D
                    counts_src = bincounts(h2_src)[:, j]
                    counts_bkg = bincounts(h2_bkg)[:, j]

                    @assert all(unique(diff(binedges(h2_bkg)[1])) .≈ Δt)
                    @assert all(unique(diff(binedges(h2_src)[1])) .≈ Δt)
                    t_bin = binedges(h2_src)[1][1:(end - 1)]

                else
                    counts_src = bincounts(h2_src[j])
                    counts_bkg = bincounts(h2_bkg[j])
                    @assert all(unique(diff(binedges(h2_bkg[j]))) .≈ Δt)
                    @assert all(unique(diff(binedges(h2_src[j]))) .≈ Δt)
                    t_bin = binedges(h2_src[j])[1:(end - 1)]

                end
                energy = "$(energies[j])-$(energies[j + 1])"
            end
            data[instr][energy] = Dict()

            # extract good time intervals
            gti_src = hdu_src[GTI_extension(observatory, CCDNR_format(curr_CCDNR_src))]
            gti_bkg = hdu_bkg[GTI_extension(observatory, CCDNR_format(curr_CCDNR_bkg))]

            # from the GTI get the Bad time intervals
            BTI_src = get_bad_time_intervals(read_header(gti_src), DataFrame(gti_src), bin_list[1], bin_list[end])
            BTI_bkg = get_bad_time_intervals(read_header(gti_bkg), DataFrame(gti_bkg), bin_list[1], bin_list[end])
            # Compute the fraction exposure for each bin
            Frac_EXP_src = get_fractional_exposure(t_bin, BTI_src, Δt)
            Frac_EXP_bkg = get_fractional_exposure(t_bin, BTI_bkg, Δt)

            if Frac_EXP_src[end] > 1.0
                Frac_EXP_src[end] = 1.0
            end
            if Frac_EXP_bkg[end] > 1.0
                Frac_EXP_bkg[end] = 1.0
            end

            println(Frac_EXP_src)
            # sanity checks
            @assert all(0.0 .<= Frac_EXP_src .<= 1.0) "Some values in Frac_EXP_src are outside the interval [0,1], $(extrema(Frac_EXP_src))"
            @assert all(0.0 .<= Frac_EXP_bkg .<= 1.0) "Some values in Frac_EXP_bkg are outside the interval [0,1], $(extrema(Frac_EXP_bkg))"

            # When fracexp is between min_Frac_EXP and 1 (not included), rescale the counts

            q_src_FEXP = min_Frac_EXP .<= Frac_EXP_src .< 1.0
            q_bkg_FEXP = min_Frac_EXP .<= Frac_EXP_bkg .< 1.0

            # arrays for the counts and variance on the counts
            u_src = float(copy(counts_src))
            u_src2 = float(copy(counts_src))
            u_bkg = float(copy(counts_bkg))
            u_bkg2 = float(copy(counts_bkg))

            # rescale the counts variance
            u_src2[q_src_FEXP] = u_src[q_src_FEXP] ./ (Frac_EXP_src[q_src_FEXP] .^ 2)
            u_bkg2[q_bkg_FEXP] = u_bkg[q_bkg_FEXP] ./ (Frac_EXP_bkg[q_bkg_FEXP] .^ 2)
            # rescale the counts

            u_src[q_src_FEXP] = u_src[q_src_FEXP] ./ Frac_EXP_src[q_src_FEXP]
            u_bkg[q_bkg_FEXP] = u_bkg[q_bkg_FEXP] ./ Frac_EXP_bkg[q_bkg_FEXP]

            # these are masks for the "good" bins
            q_src = (Frac_EXP_src .>= min_Frac_EXP)
            q_bkg = (Frac_EXP_bkg .>= min_Frac_EXP)

            # when 0<fracexp < min_Frac_EXP we can replace the data by some interpolated value and randomise
            if interpolate && observatory != :nustar
                x_interp_src = interpolate_randomise_counts(rng, t_bin, counts_src, q_src)
                x_interp_bkg = interpolate_randomise_counts(rng, t_bin, counts_bkg, q_bkg)
                @info "Interpolating and randomising small gaps number of bins filled: $(sum(.!q_src))"
            else
                u_src = allowmissing(u_src)
                u_bkg = allowmissing(u_bkg)

                u_src2 = allowmissing(u_src2)
                u_bkg2 = allowmissing(u_bkg2)

                x_interp_src, x_interp_bkg = Vector{Union{Float64, Missing}}(missing, sum(.!q_src)), Vector{Union{Float64, Missing}}(missing, sum(.!q_bkg))
            end
            # replace the bad counts
            u_src[.!q_src] = x_interp_src
            u_bkg[.!q_bkg] = x_interp_bkg
            u_src2[.!q_src] = x_interp_src
            u_bkg2[.!q_bkg] = x_interp_bkg


            data[instr][energy]["net"] = u_src - region_scale * u_bkg
            data[instr][energy]["src"] = u_src
            data[instr][energy]["bkg"] = u_bkg * region_scale

            data[instr][energy]["var_src"] = u_src2
            data[instr][energy]["var_net"] = u_src2 + region_scale^2 * u_bkg
            data[instr][energy]["var_bkg"] = u_bkg2 * region_scale^2
            data[instr]["src_interpolated"] = .!q_src
            data[instr]["bkg_interpolated"] = .!q_bkg
            data[instr]["src_fracexp"] = Frac_EXP_src
            data[instr]["bkg_fracexp"] = Frac_EXP_bkg
            data[instr]["scale"] = region_scale
            data[instr]["indexes"] = index_observations
            data[instr]["src_BTI"] = BTI_src #.-min_t_start
            data[instr]["bkg_BTI"] = BTI_bkg #.-min_t_start
        end
    end
    data["time"] = t_bin #.-min_t_start
    return data
end

function create_DimArray_from_data(data, energies, instruments)
    t = data["time"]
    T0 = t[1]
    dimt = Dim{:time}(t .- T0)
    dimen = Dim{:Erange}(get_energy_bands(energies))
    diminstr = Dim{:instrument}(Symbol.(instruments))
    dimfrac = Dim{:fracexp}([:srcfracexp, :bkgfracexp])
    dimindex = Dim{:quality}([:indexes, :interpolated])

    df_net = DimArray(reshape(mapreduce(permutedims, hcat, [mapreduce(permutedims, vcat, [data[instr][en]["net"] for instr in instruments]) for en in get_energy_bands(energies)]), (length(t), length(instruments), length(energies))), (dimt, diminstr, dimen), name = :net)
    df_varnet = DimArray(reshape(mapreduce(permutedims, hcat, [mapreduce(permutedims, vcat, [data[instr][en]["var_net"] for instr in instruments]) for en in get_energy_bands(energies)]), (length(t), length(instruments), length(energies))), (dimt, diminstr, dimen), name = :varnet)
    df_src = DimArray(reshape(mapreduce(permutedims, hcat, [mapreduce(permutedims, vcat, [data[instr][en]["src"] for instr in instruments]) for en in get_energy_bands(energies)]), (length(t), length(instruments), length(energies))), (dimt, diminstr, dimen), name = :src)
    df_bkg = DimArray(reshape(mapreduce(permutedims, hcat, [mapreduce(permutedims, vcat, [data[instr][en]["bkg"] for instr in instruments]) for en in get_energy_bands(energies)]), (length(t), length(instruments), length(energies))), (dimt, diminstr, dimen), name = :bkg)
    df_varbkg = DimArray(reshape(mapreduce(permutedims, hcat, [mapreduce(permutedims, vcat, [data[instr][en]["var_bkg"] for instr in instruments]) for en in get_energy_bands(energies)]), (length(t), length(instruments), length(energies))), (dimt, diminstr, dimen), name = :varbkg)

    df_timedata = DimArray(reshape([mapreduce(permutedims, vcat, [data[instr]["indexes"][1:(end - 1)] for instr in instruments]) mapreduce(permutedims, vcat, [data[instr]["src_interpolated"] for instr in instruments])]', (length(t), 2, length(instruments))), (dimt, dimindex, diminstr), name = :timebins)
    df_frac = DimArray(reshape([mapreduce(permutedims, vcat, [data[instr]["src_fracexp"] for instr in instruments]) mapreduce(permutedims, vcat, [data[instr]["bkg_fracexp"] for instr in instruments])]', (length(t), 2, length(instruments))), (dimt, dimfrac, diminstr), name = :frac)

    return DimStack(df_net, df_varnet, df_src, df_bkg, df_varbkg, df_frac, df_timedata), T0
end


"""
    filter_DimArray_for_overlap(dm)

    Filter the DimArray to ensure that the instruments are observing at the same time.
    This effectively truncates the
"""
function filter_DimArray_for_overlap(dm, T0, Δt_user, prefix = "")
    t = dims(dm, :time)

    file = open("$(prefix)dm.csv", "w")
    write(file, "T0:$T0\ntime,instrument,Erange,fracexp,quality,net,varnet,src,bkg,varbkg,frac,timebins\n")
    CSV.write(file, dm, append = true)
    close(file)

    # we want to find the min and max time where the instruments are all on:
    N = length(t)
    simultaneous = sum(dm.timebins[quality = At(:indexes)], dims = :instrument)
    min_index, max_index = findmax(first, simultaneous)[2][1], N - findmax(first, reverse(simultaneous))[2][1]

    dm_filtered = dm[time = t[min_index] .. t[max_index]]

    file = open("$(prefix)dm_filtered.csv", "w")
    write(file, "T0:$T0\ntime,instrument,Erange,fracexp,quality,net,varnet,src,bkg,varbkg,frac,timebins\n")
    CSV.write(file, dm_filtered, append = true)
    close(file)

    netcr = sum(dm_filtered.net, dims = :instrument) / Δt_user
    err_netcr = sqrt.(sum(dm_filtered.varnet, dims = :instrument)) / Δt_user
    bkgcr = sum(dm_filtered.bkg, dims = :instrument) / Δt_user
    err_bkgcr = sqrt.(sum(dm_filtered.varbkg, dims = :instrument)) / Δt_user

    df_err_net = DimArray(err_netcr, name = :err_net)
    df_net = DimArray(netcr, name = :net)
    df_err_bkg = DimArray(err_bkgcr, name = :err_bkg)
    df_bkg = DimArray(bkgcr, name = :bkg)

    dm_clean = DimStack(df_net, df_err_net, df_bkg, df_err_bkg)

    file = open("$(prefix)dm_clean.csv", "w")
    write(file, "T0:$T0\ntime,instrument,Erange,net,err_net,bkg,err_bkg\n")
    CSV.write(file, dm_clean, append = true)
    close(file)
    return dm_filtered, dm_clean
end
