AM_CFLAGS = $(GSL_CFLAGS) $(CFITSIO_CFLAGS) $(WCSTOOLS_CFLAGS)
AM_LDFLAGS = $(GSL_LIBS) $(CFITSIO_LIBS) $(WCSTOOLS_LIBS)

suppl = spc_driz.c spc_spc.c spc_utils.c spc_sex.c \
	spc_CD.c trace_conf.c spc_FITScards.c spc_resp.c \
	disp_conf.c aper_conf.c spc_extract.c spce_is_in.c \
	spce_sect.c spc_trace_functions.c spce_pathlength.c \
	aXe_errors.c aXe_utils.c spce_binning.c \
	spce_output.c spc_wl_calib.c spc_back.c \
	spc_flatfield.c spc_cfg.c inout_aper.c \
	spce_fitting.c spce_PET.c aper_check.c \
	drizzle_utils.c drz2pet_utils.c crossdisp_utils.c \
	spc_model.c model_utils.c spc_fluxcube.c \
	spc_optimum.c fringe_conf.c fringe_model.c \ \
	fringe_utils.c trfit_utils.c lm_eval.c lmmin.c \
	ipixcorr_utils.c inima_utils.c nicback_utils.c \
	specmodel_utils.c dirimage_model.c scaleback_utils.c

bin_PROGRAMS = aXe_SEX2GOL \
	       aXe_GOL2AF \
	       aXe_AF2PET \
	       aXe_BE \
	       aXe_PET2SPC \
	       aXe_STAMPS \
	       aXe_DRZPREP \
	       aXe_PETCONT \
	       aXe_PETFF \
	       aXe_DRZ2PET \
	       aXe_GPS \
	       aXe_FILET \
	       aXe_FRIGEN \
	       aXe_FRINGECORR \
	       aXe_TFIT \
	       aXe_INTPIXCORR \
	       aXe_PETIPC \
	       aXe_NICBACK \
	       aXe_TEST \
	       aXe_DIRIMAGE \
	       aXe_SCALEBCK

aXe_SEX2GOL_SOURCES = aXe_SEX2GOL.c $(suppl)
aXe_GOL2AF_SOURCES = aXe_GOL2AF.c $(suppl)
aXe_AF2PET_SOURCES = aXe_AF2PET.c $(suppl)
aXe_BE_SOURCES = aXe_BE.c $(suppl)
aXe_PET2SPC_SOURCES = aXe_PET2SPC.c $(suppl)
aXe_STAMPS_SOURCES = aXe_STAMPS.c $(suppl)
aXe_DRZPREP_SOURCES = aXe_DRZPREP.c $(suppl)
aXe_PETCONT_SOURCES = aXe_PETCONT.c $(suppl)
aXe_PETFF_SOURCES = aXe_PETFF.c $(suppl)
aXe_DRZ2PET_SOURCES = aXe_DRZ2PET.c $(suppl)
aXe_GPS_SOURCES = aXe_GPS.c $(suppl)
aXe_FILET_SOURCES = aXe_FILET.c $(suppl)
aXe_FRIGEN_SOURCES = aXe_FRIGEN.c $(suppl)
aXe_FRINGECORR_SOURCES = aXe_FRINGECORR.c $(suppl)
aXe_TFIT_SOURCES = aXe_TFIT.c $(suppl)
aXe_INTPIXCORR_SOURCES = aXe_INTPIXCORR.c $(suppl)
aXe_PETIPC_SOURCES = aXe_PETIPC.c $(suppl)
aXe_NICBACK_SOURCES = aXe_NICBACK.c $(suppl)
aXe_TEST_SOURCES = aXe_TEST.c $(suppl)
aXe_DIRIMAGE_SOURCES = aXe_DIRIMAGE.c $(suppl)
aXe_SCALEBCK_SOURCES = aXe_SCALEBCK.c $(suppl)
