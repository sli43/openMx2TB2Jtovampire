fftw_include=/opt/homebrew/include
fftw_lib=/opt/homebrew/lib -lfftw3

All:
	gfortran -ffree-line-length-312 -I$(fftw_include) -L$(fftw_lib) -o spw compute_dynamic.f90
