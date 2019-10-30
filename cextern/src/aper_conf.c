/*
 * Subroutines that deal with the aXe configuration file.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include "aper_conf.h"
#include "disp_conf.h"
#include "aXe_grism.h"
#include "spc_cfg.h"
#include "disp_conf.h"
#include "aXe_utils.h"

#define MAX(x,y) (((x)>(y))?(x):(y))

/**
 *  Function: get_aperture_descriptor
 *  Read a configuration file and populate an aperture_conf structure
 *  containing a description of the optical set up.
 *  This function populates the geometrical description of all the beams in
 *  all the apertures listed in the configuration file.
 *
 *  Parameters:
 *  @param actinfo  - the header structure with the column names
 *  @param filename - the complete filename to the axe configuration file
 *
 *  Returns:
 *  @return config - the configuration structure created
 */
aperture_conf *
get_aperture_descriptor (char *filename)
{
  char beam[MAXCHAR] = "\0";
  aperture_conf *config;
  int i, ix;
  gsl_vector *v;

  struct CfgStrings AperConfig[] = {
    {"INSTRUMENT", NULL},
    {"CAMERA", NULL},
    {"SCIENCE_EXT", NULL},
    {"ERRORS_EXT", NULL},
    {"DQ_EXT",NULL},
    {"DQMASK",NULL},
    {"EXPTIME",NULL},
    {"GAIN",NULL},
    {"OPTKEY1",NULL},
    {"OPTVAL1",NULL},
    {"FFNAME",NULL},
    {"REFX",NULL},
    {"REFY",NULL},
    /* aXe1.4:*/
    {"DRZRESOLA",NULL},
    {"DRZSCALE",NULL},
    {"DRZLAMB0",NULL},
    {"DRZXINI",NULL},
    {"DRZPFRAC",NULL},
    {"DRZKERNEL",NULL},
    /* aXe1.5: */
    {"PSFCOEFFS",NULL},
    {"PSFRANGE",NULL},
    {"RDNOISE",NULL},
    /* NICMOS HLA */
    {"IPIXFUNCTION",NULL},
    {"POBJSIZE",NULL},
    {"SMFACTOR",NULL},

    {NULL, NULL},
    {NULL, NULL}                /* array terminator. REQUIRED !!! */
  };

  struct CfgStrings BeamConfig[] = {
    {NULL, NULL},
    {NULL, NULL}
  };

  AperConfig[25].name = beam;
  BeamConfig[0].name = beam;


  config = malloc (sizeof (aperture_conf));
  if (config == NULL)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "Could not allocate memory for aperture configuration");
    }

  CfgRead (filename, AperConfig);

  /* Initialize optional keywords to default values*/
  sprintf(config->instrument,"None");
  sprintf(config->camera,"None");
  sprintf(config->science_ext,"None");
  sprintf(config->errors_ext,"None");
  sprintf(config->dq_ext,"None");
  sprintf(config->exptimekey,"None");
  sprintf(config->gainkey,"None");
  config->dqmask=0;
  sprintf(config->optkey1,"None");
  sprintf(config->optval1,"None");
  sprintf(config->FFname,"None");
  sprintf(config->IPIXfunc,"None");
  sprintf(config->drz_kernel,"square");
  config->refx       = 0;
  config->refy       = 0;
  config->drz_resol  = 0.0;
  config->drz_scale  = 0.0;
  config->drz_lamb0  = 0.0;
  config->drz_xstart = 15.0;
  config->drz_pfrac  = 1.0;
  config->psfcoeffs  = NULL;
  config->psfrange   = NULL;
  config->rdnoise    = 0.0;
  config->pobjsize   = -1.0;
  config->smfactor   = -1.0;

  for (ix = 0; ix < 26; ix++)
    {

      /* Name of the instrument */
      if (!strcmp (AperConfig[ix].name, "INSTRUMENT"))
        {
          if (AperConfig[ix].data == NULL)
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "INSTRUMENT tag was not read from %s\n",
                         filename);
          sprintf (config->instrument, "%s", AperConfig[ix].data);
        }

      /* Name of the camera */
      if (!strcmp (AperConfig[ix].name, "CAMERA"))
        {
          if (AperConfig[ix].data == NULL)
            aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                         "camera tag was not read from %s\n",
                         filename);
          sprintf (config->camera, "%s", AperConfig[ix].data);
        }

      /* Name of the Science data extension */
      if (!strcmp (AperConfig[ix].name, "SCIENCE_EXT"))
        {
          if (AperConfig[ix].data != NULL)
            sprintf (config->science_ext, "%s", AperConfig[ix].data);
        }

      /* Name of the optional Error FITS extension */
      if (!strcmp (AperConfig[ix].name, "ERRORS_EXT"))
        {
          if (AperConfig[ix].data != NULL)
            {
              sprintf (config->errors_ext, "%s", AperConfig[ix].data);
            }
        }

      /* Name of the optional Data Quality FITS extension */
      if (!strcmp (AperConfig[ix].name, "DQ_EXT"))
        {
          if (AperConfig[ix].data != NULL)
            {
              sprintf (config->dq_ext, "%s", AperConfig[ix].data);
            }
        }

      /* Name of the optional EXPTIME FITS keyword */
      if (!strcmp (AperConfig[ix].name, "EXPTIME"))
        {
          if (AperConfig[ix].data != NULL)
            {
              sprintf (config->exptimekey, "%s", AperConfig[ix].data);
            }
        }

      /* Red in the readout noise */
      if (!strcmp (AperConfig[ix].name, "RDNOISE"))
        {
          if (AperConfig[ix].data != NULL)
            {
              config->rdnoise = atof(AperConfig[ix].data);
            }
        }

      /* Name of the optional GAIN FITS keyword */
      if (!strcmp (AperConfig[ix].name, "GAIN"))
        {
          if (AperConfig[ix].data != NULL)
            {
              sprintf (config->gainkey, "%s", AperConfig[ix].data);
            }
        }

      /* Gett the Data Quality Mask */
      if (!strcmp (AperConfig[ix].name, "DQMASK"))
        {
          if (AperConfig[ix].data != NULL)
            {
              config->dqmask = atoi(AperConfig[ix].data);
            }
        }


      /* Optional OPTional 1st FITS keyword selection */
      if (!strcmp (AperConfig[ix].name, "OPTKEY1"))
        {
          if (AperConfig[ix].data != NULL)
            {
              sprintf (config->optkey1, "%s", AperConfig[ix].data);
            }
        }

      /* Optional OPTional 1st FITS keyword value */
      if (!strcmp (AperConfig[ix].name, "OPTVAL1"))
        {
          if (AperConfig[ix].data != NULL)
            {
              sprintf (config->optval1, "%s", AperConfig[ix].data);
            }
        }

      /* Optional FFname keyword value */
      if (!strcmp (AperConfig[ix].name, "FFNAME"))
        {
          if (AperConfig[ix].data != NULL)
            {
              sprintf (config->FFname, "%s", AperConfig[ix].data);
            }
        }

      /* Optional intra pixel correction keyword value */
      if (!strcmp (AperConfig[ix].name, "IPIXFUNCTION"))
        {
          if (AperConfig[ix].data != NULL)
            {
              sprintf (config->IPIXfunc, "%s", AperConfig[ix].data);
            }
        }

      /* Optional REFX keyword value */
      if (!strcmp (AperConfig[ix].name, "REFX"))
        {
          if (AperConfig[ix].data != NULL)
            {
              config->refx = atoi(AperConfig[ix].data);
            }
        }

      /* Optional REFY keyword value */
      if (!strcmp (AperConfig[ix].name, "REFY"))
        {
          if (AperConfig[ix].data != NULL)
            {
              config->refy = atoi(AperConfig[ix].data);
            }
        }
      if (!strcmp (AperConfig[ix].name, "DRZRESOLA"))
        {
          if (AperConfig[ix].data != NULL)
            {
              config->drz_resol = atof(AperConfig[ix].data);
            }
        }
      if (!strcmp (AperConfig[ix].name, "DRZSCALE"))
        {
          if (AperConfig[ix].data != NULL)
            {
              config->drz_scale = atof(AperConfig[ix].data);
            }
        }
      if (!strcmp (AperConfig[ix].name, "DRZLAMB0"))
        {
          if (AperConfig[ix].data != NULL)
            {
              config->drz_lamb0 = atof(AperConfig[ix].data);
            }
        }
      if (!strcmp (AperConfig[ix].name, "DRZXINI"))
        {
          if (AperConfig[ix].data != NULL)
            {
              config->drz_xstart = atof(AperConfig[ix].data);
            }
        }
      if (!strcmp (AperConfig[ix].name, "DRZPFRAC"))
        {
          if (AperConfig[ix].data != NULL)
            {
              config->drz_pfrac = atof(AperConfig[ix].data);
            }
        }
      if (!strcmp (AperConfig[ix].name, "DRZKERNEL"))
        {
          if (AperConfig[ix].data != NULL)
            {
              sprintf (config->drz_kernel, "%s", AperConfig[ix].data);
            }
        }

      if (!strcmp (AperConfig[ix].name, "PSFCOEFFS"))
        {
          if (AperConfig[ix].data != NULL)
            {
              //              sprintf (config->drz_kernel, "%s", AperConfig[ix].data);
              v = string_to_gsl_array (AperConfig[ix].data);
              config->psfcoeffs = v;
            }
        }

      if (!strcmp (AperConfig[ix].name, "PSFRANGE"))
        {
          if (AperConfig[ix].data != NULL)
            {
              v = string_to_gsl_array (AperConfig[ix].data);
              if (v->size != 2)
                aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                             "configuration file: two items for PSFRANGE, lambda_min and lambda_max must be given, not %i\n", v->size);
              config->psfrange = v;
            }
        }

      // read in the size of point-like objects
      if (!strcmp (AperConfig[ix].name, "POBJSIZE"))
        {
          if (AperConfig[ix].data != NULL)
            {
              config->pobjsize = atof(AperConfig[ix].data);
            }
        }

      // read in the adjustment for smoothed flux conversion
      if (!strcmp (AperConfig[ix].name, "SMFACTOR"))
        {
          if (AperConfig[ix].data != NULL)
            {
              config->smfactor = atof(AperConfig[ix].data);
            }
        }


    }

  // release memory
  i=0;
  while(AperConfig[i].name!=NULL)
    free(AperConfig[i++].data);

  /* Looking for up to MAX_BEAMS beams */
  config->nbeams = 0;
  for (i = 0; i < MAX_BEAMS; i++)
    {
      sprintf (beam, "BEAM%c", BEAM (i));
      CfgRead (filename, BeamConfig);
      if (BeamConfig[0].data != NULL)
        {
          v = string_to_gsl_array (BeamConfig[0].data);
          if (v->size != 2)
            {
              aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                           "tag %s does not have the proper format in file %s.\n"
                           "Should be of the form %s [integer] [integer]\n",
                           BeamConfig[0].data, filename,beam);
              //                           beam, filename,beam);
            }
          (config->beam[config->nbeams]).offset.dx0 =
            (int) gsl_vector_get (v, 0);
          (config->beam[config->nbeams]).offset.dx1 =
            (int) gsl_vector_get (v, 1);

          // release the memory
          gsl_vector_free(v);

          (config->beam[config->nbeams]).ID = i;

          // release memory
          ix=0;
          while(BeamConfig[ix].name!=NULL)
            free(BeamConfig[ix++].data);
          BeamConfig[0].data = NULL;

          (config->beam[config->nbeams]).mmag_extract = get_beam_mmag_extract (filename, i);
          (config->beam[config->nbeams]).mmag_mark = get_beam_mmag_mark (filename, i);

          sprintf (beam, "PSF_OFFSET_%c", BEAM (i));
          CfgRead (filename, BeamConfig);
          if (BeamConfig[0].data != NULL)
            {
              (config->beam[config->nbeams]).psf_offset = atof(BeamConfig[0].data);
            }
          else
            {
              (config->beam[config->nbeams]).psf_offset = 0.0;
            }

          // release memory
          ix=0;
          while(BeamConfig[ix].name!=NULL)
            free(BeamConfig[ix++].data);

          BeamConfig[0].data = NULL;

          config->nbeams = config->nbeams + 1;
        }
    }

  if (config->nbeams == 0)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
                   "No beam was found in %s. At least one beam must be found.\n",
                   filename);
    }
  return config;
}

/**
 * Function: aperture_conf_fprintf
 * Function to output the content of an aperture_conf structure.
 *
 * Parameters:
 * @param file - a pointer to an output stream or file
 * @param conf - a pointer to an existing aperture_conf structure
 */
void
aperture_conf_fprintf (FILE * file, aperture_conf * conf)
{
     int i = 0, id;

     fprintf (file, "Aperture_conf: INSTRUMENT: %s\n", conf->instrument);
     fprintf (file, "Aperture_conf: NBEAMS: %d\n", conf->nbeams);
     for (i = 0; i < conf->nbeams; i++)
       {
         id = conf->beam[i].ID;
         fprintf (file, "Aperture_conf: BEAM%c: %d %d\n", BEAM (id),
                  conf->beam[i].offset.dx0, conf->beam[i].offset.dx1);
       }
     fprintf (file, "FITS SCI extension name: %s, number: %d\n",
              conf->science_ext, conf->science_numext);
     fprintf (file, "FITS ERR extension name: %s, number: %d\n",
              conf->errors_ext, conf->errors_numext);
     fprintf (file, "FITS DQ extension name: %s, number: %d\n",
              conf->dq_ext, conf->dq_numext);
     fprintf (file, "DQ mask: %d\n",conf->dqmask);
     fprintf (file, "FF filename: %s\n",conf->FFname);
     fprintf (file, "Refx: %d\n",conf->refx);
     fprintf (file, "Refy: %d\n",conf->refy);
     fprintf (file, "DRZRESO: %10.3f\n",conf->drz_resol);
     fprintf (file, "DRZSCALE: %10.3f\n",conf->drz_scale);
     fprintf (file, "DRZLAMB0: %10.3f\n",conf->drz_lamb0);
     fprintf (file, "DRZXINI: %10.3f\n",conf->drz_xstart);
     fprintf (file, "DRZPFRAC: %10.3f\n",conf->drz_pfrac);
     fprintf (file, "DRZKERNEL: %s\n",conf->drz_kernel);


}

/**
 * Function: get_extension_numbers
 * This function assigns HDU numbers to an aperture_conf
 * structues after looking for the location of named
 * extension in a given FITS file. Optional values not
 * found are set to -1.
 *
 * Parameters:
 * @param filename - a pointer to a string containing the name
 *                   of an existing FITS file
 * @param conf     - a pointer to a populated aperture_conf structure
 * @param keyword  - a pointer to a string containing the name of
 *                   an optional keyword to read
 * @param keyval   - a pointer to a string containing the value of
 *                   the optional keyword to look for
 */
void
get_extension_numbers(char filename[], aperture_conf * conf,
                      char keyword[],char keyval[])
{
  int extver = -1;

  conf->science_numext = get_hdunum_from_hduname(filename,
                                                 conf->science_ext,keyword
                                                 ,keyval,extver);
  if (conf->science_numext > -1)
    extver = (int)get_float_from_keyword(filename, conf->science_numext, "EXTVER");

  if (strcmp(conf->errors_ext,"None"))
  {
    conf->errors_numext = get_hdunum_from_hduname(filename,
                                                  conf->errors_ext,keyword,
                                                  keyval,extver);
  } else {
    conf->errors_numext = -1;
  }
  if (strcmp(conf->dq_ext,"None")) {
    conf->dq_numext = get_hdunum_from_hduname(filename,conf->dq_ext,
                                              keyword,keyval,extver);
  } else {
    conf->dq_numext = -1;
  }
}

/**
 * Function: free_aperture_conf
 * The function releases all memory allocated
 * in a configuration structure.
 *
 * Parameters:
 * @param conf - the configuration structure
 */
void
free_aperture_conf(aperture_conf * conf)
{
  if (conf->psfcoeffs)
    gsl_vector_free(conf->psfcoeffs);
  if (conf->psfrange)
    gsl_vector_free(conf->psfrange);
  free(conf);
  conf = NULL;
}


/**
 * Function: get_psf_offset
 * The function extracts and returns the psf-offset
 * stored in a given configuration structure and
 * a given beam
 *
 * Parameters:
 * @param conf - the configuration structure
 * @param beam - the beam
 *
 * Returns:
 * @return psf_offset - the psf-offset for the beam
 */
double
get_psf_offset(aperture_conf * conf, const beam actbeam)
{
  // initialize the offset
  double psf_offset=0.0;

  int i=0;

  // go over all beams in the configuration structure
  for (i = 0; i < conf->nbeams; i++)
    {
      // check whether the current beam is the
      // right one
      if (conf->beam[i].ID == actbeam.ID)
        // return the offset
        return conf->beam[i].psf_offset;
    }

  // return the default offset
  return psf_offset;
}

/**
 * Function: get_psf_offset
 * The function extracts and returns the psf-offset
 * stored in a given configuration structure and
 * a given beam
 *
 * Parameters:
 * @param conf - the configuration structure
 * @param beam - the beam
 *
 * Returns:
 * @return psf_offset - the psf-offset for the beam
 */
double
get_max_offset(aperture_conf * conf)
{
  // initialize the offset
  double max_offset=0.0;

  int i=0;

  // go over all beams in the configuration structure
  for (i = 0; i < conf->nbeams; i++)
    {
      // check whether the current beam is the
      // right one
      max_offset = MAX(conf->beam[i].psf_offset, max_offset);
    }

  // return the default offset
  return max_offset;
}


/**
 * Function: extend_config_beams
 * The function extends the size of the beams in a
 * configuration structure.
 *
 * Parameters:
 * @param conf - the configuration structure
 */
void
extend_config_beams(aperture_conf *conf)
{
  int index;

  // go over all beams in the configuration structure
  for (index=0; index < conf->nbeams; index++)
    {
      // lower the left border
      conf->beam[index].offset.dx0 -= BCK_BEAM_EXT;

      // enhance the right border
      conf->beam[index].offset.dx1 += BCK_BEAM_EXT;
    }
}
