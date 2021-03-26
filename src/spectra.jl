using FITSIO
using MonteCarloMeasurements

"""
    read_spectrum_unknown_source_prof(filepath)

Read one of the spectra given by the professor - source unknown. The FITS file is assumed 
to have a header with the following keys:

- 'CRVAL1': initial wavelength
- 'CDELT1': lambda wavelength

such that to each pixel with index `i` is assigned the wavelength ``CRVAL1 + CDELT1*i``.

# Arguments
- `filepath`: the relative or absolute filepath of the spectrum file.
"""
function read_spectrum_unknown_source_prof(filepath)
    f = FITS(filepath)
    hdu = f[1]

    header = read_header(hdu)
    flux = read(hdu)

    Δx = header["CDELT1"]
    λᵢ = header["CRVAL1"]

    λ_space = collect(λᵢ .+ (1:length(flux))*Δx)

    # TODO: Do I need a struct here?
    (λ_space, flux, Δx)
end

"""
    read_spectrum_ambre_grid(filepath)

Read a wavelength or flux file from the AMBRE GRID project.

# Arguments
- `filepath`: the relative or absolute filepath of the spectrum file.
"""
function read_spectrum_ambre_grid(filepath)
    f = FITS(filepath)
    hdu = f[1]
    read(hdu)
end

"""
    read_spectrum_sisma(filepath)

Read a SISMA/SpaceINN database spectrum.
"""
function read_spectrum_sisma(filepath)
    f = FITS(filepath)
    hdu = f[2]

    λ_space = read(hdu, "WAVE")
    flux = read(hdu, "FLUX")

    header = read_header(hdu)

    # NOTE: there should be no type instability because the numerical values
    # are all 64bit floats
    vrad = typeof(header["VRAD"]) == String ? missing : header["VRAD"]
    vrad_error = typeof(header["VRAD_ERR"]) == String ? missing : header["VRAD_ERR"]
    Teff = header["TEFF"]
    logg = header["LOGG"]
    FeH = header["FEH"]
    vsini = typeof(header["VSINI"]) == String ? missing : header["VSINI"]
    vsini_error = typeof(header["VSIN_ERR"]) == String ? missing : header["VSIN_ERR"]
    R = header["SPEC_RP"]

    opt_data = Dict(
        ("VRAD" => vrad ± vrad_error,
         "TEFF" => Teff,
         "LOGG" => logg,
         "FEH" => FeH,
         "VSINI" => vsini ± vsini_error,
         "R" => R)
    )

    (λ_space, flux, opt_data)
end

# TODO: add return values to docstrings