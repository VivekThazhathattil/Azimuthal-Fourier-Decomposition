#ifndef _H5_WRAPPER_H
#define _H5_WRAPPER_H

	#include <hdf5.h>	
	#include <stdlib.h>
	#include <complex.h>
	#include <math.h>

	typedef struct IN_DSET_s{
		int nt; 			/* number of time snapshots */
		int nth; 			/* number of azimuthal planes */
		int nrnz; 			/* total number of points in RZ plane */
		double** dat; /* data in the form (nt, nth*nr*nz) */
	} in_dset_t;

	typedef struct COMPLEX_DTYPE_s{
		double re;	 /* real part */
		double im; 	 /* imag part */
	} complex_t;

	typedef struct AZ_FR_s{
		double _Complex *ft;
		double *ift;
	} az_ft_res_t;

	typedef struct H5_SAVE_s{
		double _Complex** ft;
		double** ift;
		
	}	h5_save_t;

	typedef double _Complex DC_t;

	struct IN_DSET_s* fetch_h5_data(const char*, const char*);
	struct IN_DSET_s* arrange_data(double*, hsize_t, hsize_t, hsize_t);
	struct H5_SAVE_s* initialize_h5_save(int);
	void store_h5_data(struct AZ_FR_s*, const char*, int, int, int);	
	void create_new_dset(const char*, const char*, void*, char, hsize_t*);
	void free_1d_data(double*);
	void free_in_dset(struct IN_DSET_s*);
	void free_h5_save(struct H5_SAVE_s*);
	void save_h5(struct H5_SAVE_s*, int, int, int);
	void full_h5_save(char*, struct H5_SAVE_s*, char, hsize_t*, int);

#endif

