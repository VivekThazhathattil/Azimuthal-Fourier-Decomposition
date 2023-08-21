#include <stdio.h>
#include <omp.h>

#include "h5_wrapper.h"
#include "fourier_transform.h"

/*-------------------------------------------------------------*/
int main(int argc, char* argv[])
/*-------------------------------------------------------------*/
{
	double beg_time, end_time;
	beg_time = omp_get_wtime();

	const char* in_file_name = "in_data/ur_104p.h5";
	const char* dset_path = "/ds1";
	const char* out_file_prefix = "out_data/ur_104p_nmall/";
	const int num_recon_modes = 103;
	char out_name_buf[256];
	int i;
	int num_threads = 24;

	/* set number of threads that OpenMP will use */
	omp_set_num_threads(num_threads);

	printf("Fetching data...\n");
	in_dset_t* dset = fetch_h5_data(in_file_name, dset_path);

	int j;

	int num_modes = get_num_modes(dset->nth);
	h5_save_t* save_data = initialize_h5_save(dset->nt);

	#pragma omp parallel for
	for(i = 0; i < dset->nt; ++i){
		double loop_beg_time, loop_end_time;
		loop_beg_time = omp_get_wtime();	

		printf("Compute Fourier decomp for tshot: %d...\n",
			 i + 1);

		az_ft_res_t* res = get_azim_fourier_decomp(dset, 
			num_recon_modes, i);		
		
		loop_end_time = omp_get_wtime();	
		printf("Time taken: %lf seconds\n\n", 
			loop_end_time - loop_beg_time);

		save_data->ft[i] = res->ft;
		save_data->ift[i] = res->ift;

		free(res);
	}

	save_h5(save_data, dset->nt, num_modes, dset->nrnz, 
		out_file_prefix);

	free_h5_save(save_data, dset->nt);

	end_time = omp_get_wtime();
	printf("Total time taken: %lf seconds\n", 
		end_time - beg_time);

	return 0;
}
