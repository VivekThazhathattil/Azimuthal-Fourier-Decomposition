#include "h5_wrapper.h"

/*-------------------------------------------------------------*/
in_dset_t* fetch_h5_data(const char* in_file, const char* dpath)
/*-------------------------------------------------------------*/
{
	hid_t file_id, dset_id, dataspace_id;	
	herr_t status;
	hsize_t dims[3]; // in this case our dspace has dim Npts x Nobs

	file_id = H5Fopen(in_file, H5F_ACC_RDONLY, H5P_DEFAULT);
	dset_id = H5Dopen2(file_id, dpath, H5P_DEFAULT);
	dataspace_id = H5Dget_space(dset_id);
	H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

	double* dset = 
		(double*) malloc(dims[0] * dims[1] * dims[2] * sizeof(double));

	status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
		H5P_DEFAULT, dset);
	
	/* dims = (21534 x 52 x 515)  = (nrn x nth x nt) */
	in_dset_t* res = arrange_data(dset, dims[2], dims[1], dims[0]);		

	H5Dclose(dset_id);
	H5Sclose(dataspace_id);
	H5Fclose(file_id);

	free_1d_data(dset);

	return res;
}

/*-------------------------------------------------------------*/
in_dset_t* arrange_data(double* dat, hsize_t nt, hsize_t nth, 
	hsize_t nrnz)
/*-------------------------------------------------------------*/
{
	int i, j, k;
	in_dset_t* mat = (in_dset_t*) malloc(sizeof(in_dset_t));
	mat->nt = nt;
	mat->nth = nth;
	mat->nrnz = nrnz;
	mat->dat = (double**) malloc(sizeof(double*) * nt);

	for(i = 0; i < nt; ++i){
		mat->dat[i] = (double*) malloc(sizeof(double) * nth * nrnz);
	}

	for(i = 0; i < nrnz; ++i){
		for(j = 0; j < nth; ++j){
			for(k = 0; k < nt; ++k){
				mat->dat[k][nth * i + j] = 
					dat[k + (nt * j) + (nt * nth * i)];	
			}
		}
	}

	return mat;
}

/*-------------------------------------------------------------*/
void create_new_dset(const char* out_file_name , 
	const char* dpath, void* data, char dtype, hsize_t* dims)
/*-------------------------------------------------------------*/
{
	hid_t file_id, dset_id, type_id, dspace_id;

	H5E_BEGIN_TRY{
		file_id = H5Fcreate(out_file_name, H5F_ACC_EXCL, H5P_DEFAULT,
			H5P_DEFAULT);
	}H5E_END_TRY
	if(file_id == H5I_INVALID_HID){
		file_id = H5Fopen(out_file_name, H5F_ACC_RDWR, H5P_DEFAULT);
	}

	if(dtype == 'd'){
		type_id = H5T_NATIVE_DOUBLE;
		data = (double*) data;
	}

	else if(dtype == 'c'){
		hid_t complex_id = H5Tcreate(H5T_COMPOUND, sizeof(complex_t));
		H5Tinsert(complex_id, "real", HOFFSET(complex_t, re),
			H5T_NATIVE_DOUBLE);
		H5Tinsert(complex_id, "imag", HOFFSET(complex_t, im),
			H5T_NATIVE_DOUBLE);
		type_id = complex_id;
		data = (double _Complex*) data;

	}

	dspace_id = H5Screate_simple(1, dims, NULL);

	dset_id = H5Dcreate2(
							file_id, 			/* location identifier for file */
							dpath,				/* name of dataset to create */
						  type_id,	 		/* datatype identifier */
						  dspace_id,	 	/* dataspace identifier */
						  H5P_DEFAULT,	/* link creation property list identifer */
						  H5P_DEFAULT,	/* dataset creation property list identifer */
						  H5P_DEFAULT		/* dataset access property list identifer */
					);	
	H5Dwrite(dset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	H5Dclose(dset_id);
	H5Sclose(dspace_id);
	H5Fclose(file_id);
}

/*-------------------------------------------------------------*/
void store_h5_data(az_ft_res_t* res, const char* out_name_buf, 
	int tshot, int nrnz, int nth)
/*-------------------------------------------------------------*/
{
	hsize_t dims[1];
	char dpath_str[256];
	char out_name_new[256];

	if(nth % 2){
		dims[0] = nrnz * nth;
	}
	else{
		dims[0] = nrnz * (nth - 1);
	}
	sprintf(out_name_new, "ft_%s", out_name_buf);
	sprintf(dpath_str, "/%d", tshot);
	create_new_dset(out_name_new, dpath_str, res->ft, 'c', dims);

	sprintf(out_name_new, "ift_%s", out_name_buf);
	sprintf(dpath_str, "/%d", tshot);
	create_new_dset(out_name_new, dpath_str, res->ift, 'd', dims);
}

/*-------------------------------------------------------------*/
void full_h5_save(char* out_file_name , h5_save_t* sav, 
	char dtype, hsize_t* dims, int nt)
/*-------------------------------------------------------------*/
{
	hid_t file_id, dset_id, type_id, dspace_id;
	int i;

	H5E_BEGIN_TRY{
		file_id = H5Fcreate(out_file_name, H5F_ACC_EXCL, H5P_DEFAULT,
			H5P_DEFAULT);
	}H5E_END_TRY
	if(file_id == H5I_INVALID_HID){
		file_id = H5Fopen(out_file_name, H5F_ACC_RDWR, H5P_DEFAULT);
	}

	if(dtype == 'd'){
		type_id = H5T_NATIVE_DOUBLE;
	}

	else if(dtype == 'c'){
		hid_t complex_id = H5Tcreate(H5T_COMPOUND, sizeof(complex_t));
		H5Tinsert(complex_id, "real", HOFFSET(complex_t, re),
			H5T_NATIVE_DOUBLE);
		H5Tinsert(complex_id, "imag", HOFFSET(complex_t, im),
			H5T_NATIVE_DOUBLE);
		type_id = complex_id;
	}

	dspace_id = H5Screate_simple(1, dims, NULL);

	for(i = 0; i < nt; ++i){
		char dpath[128];	
		sprintf(dpath, "/%d", i);
		dset_id = H5Dcreate2(
								file_id, 			/* location identifier for file */
								dpath,				/* name of dataset to create */
							  type_id,	 		/* datatype identifier */
							  dspace_id,	 	/* dataspace identifier */
							  H5P_DEFAULT,	/* link creation property list identifer */
							  H5P_DEFAULT,	/* dataset creation property list identifer */
							  H5P_DEFAULT		/* dataset access property list identifer */
						);	

		if(dtype == 'c'){
			H5Dwrite(dset_id, type_id, H5S_ALL, H5S_ALL, 
				H5P_DEFAULT, sav->ft[i]);
		}
		else{
			H5Dwrite(dset_id, type_id, H5S_ALL, H5S_ALL, 
				H5P_DEFAULT, sav->ift[i]);
		}
	}

	H5Dclose(dset_id);
	H5Sclose(dspace_id);
	H5Fclose(file_id);
}

/*-------------------------------------------------------------*/
void save_h5(h5_save_t* sav, int nt, int nmodes, int nrnz,
	const char* out_file_prefix)
/*-------------------------------------------------------------*/
{
	hsize_t dims[1];
	char out_name[256];

	dims[0] = nrnz * nmodes;

	sprintf(out_name, "%sout_ft.h5", out_file_prefix);
	full_h5_save(out_name, sav, 'c', dims, nt);

	sprintf(out_name, "%sout_ift.h5", out_file_prefix);
	full_h5_save(out_name, sav, 'd', dims, nt);
}

/*-------------------------------------------------------------*/
h5_save_t* initialize_h5_save(int nt)
/*-------------------------------------------------------------*/
{
	h5_save_t* save_data = (h5_save_t*) malloc(sizeof(h5_save_t));
	save_data->ft = (DC_t**) malloc(sizeof(DC_t*) * nt);
	save_data->ift = (double**) malloc(sizeof(double*) * nt);
	return save_data;
}

/*-------------------------------------------------------------*/
void free_1d_data(double* data)
/*-------------------------------------------------------------*/
{
	free(data);
}


/*-------------------------------------------------------------*/
void free_in_dset(in_dset_t* mat)
/*-------------------------------------------------------------*/
{
	int i;
	for(i = 0; i < mat->nt; ++i)
		free(mat->dat[i]);
	free(mat);	
}

/*-------------------------------------------------------------*/
void free_h5_save(struct H5_SAVE_s* save_data, int nt)
/*-------------------------------------------------------------*/
{
	int i;

	for(i = 0; i < nt; ++i){
		free(save_data->ft[i]);	
	}
	free(save_data->ft);

	for(i = 0; i < nt; ++i){
		free(save_data->ift[i]);
	}
	free(save_data->ift);

	free(save_data);
}
