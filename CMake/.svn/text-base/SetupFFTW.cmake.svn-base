#
# This macro is for setting up a sub-project to use FFTW (float or double version)
#

macro(SetupFFTW)
  mark_as_advanced(USE_FFTWD)
  option(USE_FFTWD "Use double precision FFTW if found" ON)
  mark_as_advanced(USE_FFTWF)
  option(USE_FFTWF "Use single precision FFTW if found" ON)

  find_package(FFTW)

  # Set the FFTW_LIB variable to the version we found.
  # Prefer the double precision version if both found.

  set(FFTW_FOUND NOTFOUND)

  if(FFTWF_FOUND)
    set(FFTW_LIB ${FFTWF_LIB})
    set(FFTW_FOUND TRUE)
  endif(FFTWF_FOUND)

  if(FFTWD_FOUND)
    set(FFTW_LIB ${FFTWD_LIB})
    set(FFTW_FOUND TRUE)
  endif(FFTWD_FOUND)

  MESSAGE(FFTW_LIB: ${FFTW_LIB})

endmacro(SetupFFTW)
