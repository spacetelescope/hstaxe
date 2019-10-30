#ifndef _SPC_FLUXCUBE_H
#define _SPC_FLUXCUBE_H


typedef struct
{
  double wavelength;    // the wavelength of the flux image  
  gsl_matrix *flux;     // the pixel values in the flux image
}
flux_image;

typedef struct
{
  int xoffs;                     // x-pixel offset between flt and fluxcube image
  int yoffs;                     // y-pixel offset between flt and fluxcube image
  int n_fimage;                  // the number of flux images
  gsl_matrix_int  *segmentation; // the segmentation image
  flux_image     **fluxims;      // the flux images
  gsl_vector_int  *fimage_order; // vector with the indices of the fluximages
}
flux_cube;

extern flux_cube *
load_fluxcube(const char fcube_file[]);

extern flux_cube *
alloc_fluxcube(const int nflux);

extern int 
load_offsets(const char fcube_file[], flux_cube *fcube);

extern gsl_matrix_int *
load_segmentation(const char fcube_file[]);

extern flux_image *
load_fluximage(const char fcube_file[], int hdunum);

extern gsl_vector_int *
order_fluxims(flux_cube *fcube);

extern void
free_fluxcube(flux_cube *fcube);

extern void
free_fluximage(flux_image *fimage);

extern d_point 
flt_to_fcube_trans(const flux_cube *fcube,const  d_point point);

extern d_point 
fcube_to_flt_trans(const flux_cube *fcube,const  d_point point);

extern dirobject **
fluxcube_to_dirlist(const flux_cube *fcube, object  **oblist);

extern dirobject *
dirobject_from_segpoint(const px_point point, const object  *obj);

extern void
update_dirobject(dirobject *actdir, const px_point fcube_point);

extern void
fill_xy_offsets(dirobject **dirlist, char CONF_file[]);

extern void 
fill_fluxvalues(const flux_cube *fcube, const px_point point,
		dirobject *actdir, const int inter_type);

extern void
print_fluxcube(const flux_cube *fcube);

#endif
