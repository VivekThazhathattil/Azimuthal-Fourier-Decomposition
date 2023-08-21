#ifndef _FOURIER_TRANSFORM_H
#define _FOURIER_TRANSFORM_H

	#include <stdio.h>
	#include <complex.h>
	#include <math.h>
	#include <assert.h>

	#include "h5_wrapper.h"
	
	az_ft_res_t* get_azim_fourier_decomp(in_dset_t*, int, int);
	DC_t trap_integrate(in_dset_t*, int, int, int);
	DC_t* get_fourier_coeffs(in_dset_t*, int);
	int get_num_modes(int);
	double* get_inv_fourier_trans(DC_t*, in_dset_t*, int);
	void free_az_ft_res(az_ft_res_t*);

#endif
