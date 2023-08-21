#include "h5_wrapper.h"
#include "fourier_transform.h"

/*-------------------------------------------------------------*/
int get_num_modes(int nth)
/*-------------------------------------------------------------*/
{
	if(nth % 2){
		return nth;
	}
	return nth - 1;
}

/*-------------------------------------------------------------*/
az_ft_res_t* get_azim_fourier_decomp(in_dset_t* dset,
	int num_modes_for_recons, int tidx)
/*-------------------------------------------------------------*/
{
	az_ft_res_t* res = (az_ft_res_t*) malloc(sizeof(az_ft_res_t));
	res->ft = get_fourier_coeffs(dset, tidx);
	res->ift = get_inv_fourier_trans(res->ft, dset, 
		num_modes_for_recons);

	return res;
}

/*-------------------------------------------------------------*/
DC_t trap_integrate(in_dset_t* df, int tidx, int pidx, int m)
/*-------------------------------------------------------------*/
{

	DC_t g(int mode, int ndivs, int div_bound_idx, double val){
		return cexp((-2 * M_PI * I * mode * div_bound_idx)/ndivs) * val;
	}

	int ii;
	int ndiv = df->nth;
	double* d = df->dat[tidx];
	double dx = (2 * M_PI) / ndiv;
	DC_t res = 0.0;

	for(ii = 0; ii < ndiv - 1; ++ii){
		res += (g(m, ndiv, ii, d[pidx + ii]) + 
						g(m, ndiv, ii, d[pidx + ii])) * (dx/2);
	}

	res += (g(m, ndiv, ndiv - 1, d[pidx + ndiv - 1]) +
					g(m, ndiv, ndiv, d[pidx])) * (dx/2);

	return res;
}

/*-------------------------------------------------------------*/
DC_t* get_fourier_coeffs(in_dset_t* dset, int tidx)
/*-------------------------------------------------------------*/
{
	printf("Computing Fourier coeffs...\n");
	int num_modes, ii, jj;

	int nrnz = dset->nrnz;
	num_modes = get_num_modes(dset->nth);

	DC_t* ft = (DC_t*) malloc(sizeof(DC_t) * num_modes * nrnz);
	for(ii = 0; ii < nrnz; ++ii){
		for(jj = 0; jj < num_modes; ++jj){
			/* We need our modes to be [-num_modes/2, ..., num_modes/2]*/
			int m = jj - num_modes/2;
			ft[ii*num_modes + jj] = trap_integrate(dset, tidx, 
				ii*dset->nth, m);
		}
	}
	return ft;
}

/*-------------------------------------------------------------*/
double* get_inv_fourier_trans(DC_t* ft, in_dset_t* dset,
	 int num_recon_modes)
/*-------------------------------------------------------------*/
{
	printf("Computing reconstructed field...\n");
	assert(num_recon_modes % 2);
	assert(num_recon_modes <= dset->nth);

	int nrnz = dset->nrnz;
	int ii, jj, kk;
	int total_num_modes = get_num_modes(dset->nth);

	double* ift =
		 (double*) calloc(nrnz * total_num_modes, sizeof(double));

	for(ii = 0; ii < nrnz; ++ii){
		for(jj = 0; jj < total_num_modes; ++jj){
			for(kk = 0; kk < num_recon_modes; ++kk){
				int m = kk - num_recon_modes/2;
				int pidx = (ii * total_num_modes) +
					(total_num_modes/2) - (num_recon_modes/2) + jj;
				ift[ii * total_num_modes + jj] += creal(ft[pidx + kk] * 
					cexp(I * M_PI * m * jj/dset->nth));
			}
		}
	}

	return ift;
}

/*-------------------------------------------------------------*/
void free_az_ft_res(az_ft_res_t* res)
/*-------------------------------------------------------------*/
{
	free(res->ft);
	free(res->ift);
	free(res);
}
