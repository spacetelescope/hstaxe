/**
 */
#include        "inout_aper.h"

/**
 * Function: object_list_to_file
 * Dump an object apperture list into a file
 *
 * Parameters:
 * @param oblist           - a point to an object list
 * @param filename         - name of the file to write output
 * @param leaveout_ignored - if ignore=1 beams are left out of the aper file
 *
 * Returns:
 * @return num - number of beams written to the aperture file
 */
int
object_list_to_file (object * const *oblist, char *filename,
                     int leaveout_ignored)
{
  FILE *file;
  object *const *obp;
  const beam *b;
  char refpixel[MAXCHAR];
  char corners[MAXCHAR];
  char curve[MAXCHAR];
  char width[MAXCHAR];
  char orient[MAXCHAR];
  //char bwindow[MAXCHAR];
  char awidth[MAXCHAR];
  char bwidth[MAXCHAR];
  char aorient[MAXCHAR];
  char flux[MAXCHAR];
  char slitgeom[MAXCHAR];
  char modspec[MAXCHAR];
  char modimage[MAXCHAR];
  char ignore[MAXCHAR];
  char tmps[MAXCHAR];
  int oid, i, ii, num = 0;
  char label[80];
  double *cf;
  int wrote_aper;
  //int j = 0;

  // open the output file;
  // give error in case of problems
  file = fopen (filename, "w");
  if (file==NULL)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Could not open %s", filename);
    }

  // go over the apertures
  for (obp = oblist; *obp; obp++)
    {
      oid = (*obp)->ID;
      wrote_aper = 0;
      for (i = 0; i < (*obp)->nbeams; i++)
        {
          b = &((*obp)->beams[i]);
          if ((leaveout_ignored) && (b->ignore == 1))
            continue;

          if (!wrote_aper)
            {
              fprintf (file, "APERTURE %i\n", oid);
              wrote_aper = 1;
            }

          sprintf (label, "Aperture %i, BEAM %c", oid, BEAM (i));

          fprintf (file, "  BEAM %c\n", BEAM (i));

          /* Output REFPIXEL info */
          sprintf (refpixel,
              "     REFPIXEL%i%c %.3f %.3f"
              " \t\t\t\t;%s Reference Pixel.", oid, BEAM (i),
              b->refpoint.x, b->refpoint.y, label);

    	  /* Output CORNERS aperture info */
          sprintf (corners,
              "     CORNERS%i%c  %i %i %i %i %i %i %i %i"
              " \t\t;%s Aperture Coordinates.", oid, BEAM (i),
              b->corners[0].x, b->corners[0].y, b->corners[1].x,
              b->corners[1].y, b->corners[2].x, b->corners[2].y,
              b->corners[3].x, b->corners[3].y, label);

    	  /* Output CURVE trace description info */
          sprintf (curve,"     CURVE%i%c    %i", oid, BEAM (i), b->spec_trace->type);
          cf = b->spec_trace->data;

          for (ii=0;ii<b->spec_trace->type + 1;ii++)
            {
              sprintf(tmps," %.3e",cf[ii+1]);
              strcat(curve,tmps);
            }
          sprintf(tmps,"\t\t\t;%s Trace description.",label);
          strcat(curve,tmps);

          /* Output WIDTH width of object info */
          sprintf (width,
              "     WIDTH%i%c    %.3f \t\t\t\t\t;%s Object Half Width.",
              oid, BEAM (i), b->width, label);

    	  /* Output ORIENT orientation angle of object info */
          /* Angle reference frame is converted from aXe's to SeXtractor */
          sprintf (orient,
              "     ORIENT%i%c   %.3f \t\t\t\t\t;%s Object Orientation.",
              oid, BEAM (i), ((b->orient) / M_PI * 180. - 180.), label);


          if (b->slitgeom[0] > -1.0)
            {
              /* Output CORNERS aperture info */
              sprintf (slitgeom,
                  "     SLITGEOM%i%c %.2f %.2f %.2f %.2f \t\t; %s Object Orientation.", oid, BEAM (i),
                  b->slitgeom[0],
                  b->slitgeom[1] / M_PI * 180. - 180.,
                  b->slitgeom[2],
                  b->slitgeom[3],
                  label);
            }

          if (b->modspec >-1)
            {
              // outputmodel spectrum value
              sprintf (modspec,
                  "     MODSPEC%i%c  %3i \t\t\t\t;%s, Index for model spectrum.",
                  oid, BEAM (i), b->modspec, label);
    	    }
          if (b->modimage >-1)
            {
              // output model image value
              sprintf (modimage,
                  "     MODIMAGE%i%c  %3i \t\t\t\t;%s, Index for direct emission.",
                  oid, BEAM (i), b->modimage, label);
    	    }

          if (b->awidth > 0.0)
            {
              /* Output AWIDTH width of object info */
              sprintf (awidth,
                  "     AWIDTH%i%c   %.3f \t\t\t\t\t\t;%s Maximum Object Half Width.(SExtractor A_IMAGE)",
                  oid, BEAM (i), b->awidth, label);

              /* Output BWIDTH width of object info */
              sprintf (bwidth,
                  "     BWIDTH%i%c   %.3f \t\t\t\t\t\t;%s Minimum Object Half Width.(SExtractir B_IMAGE)",
                  oid, BEAM (i), b->bwidth, label);

              /* Output AORIENT angle of object info */
              sprintf (aorient,
                  "     AORIENT%i%c  %.3f \t\t\t\t\t\t;%s Orientation of Object.(SExtractor THETA_IMAGE)",
                  oid, BEAM (i), b->aorient, label);

              // print the identifier
              sprintf (flux,"     FLUX%i%c    ", oid, BEAM (i));

              // go over all flux vector values
              for (ii=0;ii<(int)b->flux->size;ii++)
                {
                  // print and append the current value
                  sprintf(tmps," %.5e", gsl_vector_get(b->flux, ii));
                  strcat(flux,tmps);
                }
              // print an append the comment
              sprintf(tmps,"\t\t;%s Flux values.",label);
              strcat(flux,tmps);
            }

          /* Output IGNORE flag */
          sprintf (ignore,
              "     IGNORE%i%c   %1d  \t\t\t\t\t;%s ignore this beam.",
              oid, BEAM (i), b->ignore, label);

          // write the beam components to the file
          fprintf (file, "%s\n", refpixel);
          fprintf (file, "%s\n", corners);
          fprintf (file, "%s\n", curve);
          fprintf (file, "%s\n", width);
          fprintf (file, "%s\n", orient);
          if (b->slitgeom[0] > -1.0)
            fprintf (file, "%s\n", slitgeom);
          if (b->modspec > -1)
            fprintf (file, "%s\n", modspec);
          if (b->modimage > -1)
              fprintf (file, "%s\n", modimage);
          if (b->awidth > 0.0)
            {
              fprintf (file, "%s\n", awidth);
              fprintf (file, "%s\n", bwidth);
              fprintf (file, "%s\n", aorient);
              fprintf (file, "%s\n", flux);
            }
          fprintf (file, "%s\n", ignore);

          // mark the beam end
          fprintf (file, "  BEAM END\n");

          // enhance the beam counter
          num++;
        }
      if (wrote_aper)
        fprintf (file, "APERTURE END\n\n");
      }
  fclose (file);
  return num;
}


/**
 * Function: aper_file_apernum
 * A helper function which returns the number of time the keyword
 * APERTURE appears in an Aperture File and returns that number.
 *
 * Parameters:
 * @param name of an Aperture File
 *
 * Returns:
 * @return number of time APERTURE appears in the file
 */
int
aper_file_apernum (char *filename)
{
  int n = 0;
  char Buffer[BUFFERSIZE] = "\0";
  char *WorkPtr           = NULL;
  char *CfgName           = NULL;
  char *CfgData           = NULL;
  FILE *CfgFile;

  CfgFile = fopen (filename, "r");
  if (NULL == CfgFile)
    aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Could not open %s", filename);


  while (NULL != fgets (Buffer, BUFFERSIZE, CfgFile))
    {
      /* clip off optional comment tail indicated by a semi-colon
       */
      if (NULL != (WorkPtr = strchr (Buffer, ';')))
        *WorkPtr = '\0';
      else
        WorkPtr = Buffer + strlen (Buffer);

      /* clip off trailing and leading white space
       */
      WorkPtr--;
      while (isspace ((int) *WorkPtr) && WorkPtr >= Buffer)
        *WorkPtr-- = '\0';
      WorkPtr = Buffer;
      while (isspace ((int) *WorkPtr))
        WorkPtr++;
      if (0 == strlen (WorkPtr))
        continue;

      CfgName = strtok (WorkPtr, " =");
      if (NULL != CfgName)
        {
          /* Condition the name (lower case required),
             and strip leading white and a 'late' = from data part.
           */
          CfgData = strtok (NULL, "");
          if (CfgData != NULL)
            {
              while (isspace ((int) *CfgData))
                CfgData++;
              if ('=' == *CfgData)
                CfgData++;
              while (isspace ((int) *CfgData))
                CfgData++;
            }
        }
      if (!strcmp ("APERTURE", CfgName)  && strcmp ("END", CfgData) )
        n++;
    }
  fclose(CfgFile);
  return n;
}

/**
 * Function: aper_file_aperlist
 * This function returns a gsl_vector containing the list of apertures
 * found in the given aperture file.
 *
 * Parameters:
 * @param filename - the name of the aperture file
 *
 * Returns:
 * @return apers - a gsl_vector_int * containing the list of aperture IDs
 */
gsl_vector_int * aper_file_aperlist(char *filename)
  {
    int n = 0;
    char Buffer[BUFFERSIZE] = "\0";
    char *WorkPtr= NULL;
    char *CfgName= NULL;
    char *CfgData= NULL;
    FILE *CfgFile= NULL;
    int napers;
    gsl_vector_int *apers;

    napers = aper_file_apernum(filename);
    if (!(napers>0))
      return NULL;

    apers = gsl_vector_int_alloc(napers);

    CfgFile = fopen(filename, "r");
    if (NULL == CfgFile)
      {
        aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
        "Could not open aperture file %s\n", filename);
      }

    while (NULL != fgets (Buffer, BUFFERSIZE, CfgFile))
      {
        /* clip off optional comment tail indicated by a semi-colon
         */
        if (NULL != (WorkPtr = strchr (Buffer, ';')))
          *WorkPtr = '\0';
        else
          WorkPtr = Buffer + strlen(Buffer);

        /* clip off trailing and leading white space
         */
        WorkPtr--;
        while (isspace ((int) *WorkPtr) && WorkPtr >= Buffer)
          *WorkPtr-- = '\0';
        WorkPtr = Buffer;
        while (isspace ((int) *WorkPtr))
          WorkPtr++;
        if (0 == strlen(WorkPtr))
          continue;

        CfgName = strtok(WorkPtr, " =");
        CfgData = NULL;
        if (NULL != CfgName)
          {
            /* Condition the name (lower case required),
             and strip leading white and a 'late' = from data part.
             */
            //strlwr( CfgName );
            CfgData = strtok(NULL, "");
            if (CfgData != NULL)
              {
                while (isspace ((int) *CfgData))
                  CfgData++;
                if ('=' == *CfgData)
                  CfgData++;
                while (isspace ((int) *CfgData))
                  CfgData++;
              }
          }
        if (!strcmp("APERTURE", CfgName) && strcmp("END", CfgData))
          {
            gsl_vector_int_set(apers, n, atoi(CfgData));
            n++;
            if (n > napers)
              {
                aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                "File format error in aperture file %s",
                filename);
              }
          }
      }
    fclose(CfgFile);
    return apers;
  }

/**
 * Function: nbeams_from_aper_file
 * Read an Aperture File and returns the number of beams
 * found for a given aperture number (up to MAX_BEAMS)
 *
 * Parameters:
 * @param filename - name of Aperture File
 * @param num      - number of apertures found in Aperture File
 *
 * Returns:
 * @returns nbeam - the number of beams for the aperture
 */
int
nbeams_from_aper_file (char *filename, int num)
  {
    int i, nbeam = 0;
    char refpixel[MAXCHAR];

    for (i = 0; i < MAX_BEAMS; i++)
      {
        struct CfgStrings AperData[] =
          {
            { refpixel, NULL },
            { NULL, NULL } /* array terminator. REQUIRED !!! */
          };
        sprintf(refpixel, "REFPIXEL%d%c", num, BEAM (i));
        CfgRead(filename, AperData);

        if (!strcmp(AperData[0].name, refpixel))
          if (AperData[0].data != NULL)
            nbeam++;
      }
    return nbeam;
  }

/**
 * Function: nbeams_from_char_array2
 * Read an Aperture File and returns the number of beams
 * found for a given aperture number (up to MAX_BEAMS)
 *
 * Parameters:
 * @param filename - name of Aperture File
 * @param num      - number of apertures found in Aperture File
 *
 * Returns:
 * @return nbeam - the number of beams for the aperture
 *
 */
int
nbeams_from_char_array2 (char **apers, int num)
  {
    int i, j, nbeam = 0;
    char refpixel[MAXCHAR];

    for (i = 0; i < MAX_BEAMS; i++)
      {
        struct CfgStrings AperData[] = {
              {refpixel, NULL},
              {NULL, NULL}    /* array terminator. REQUIRED !!! */
        };
        sprintf (refpixel, "REFPIXEL%d%c", num, BEAM (i));
        CfgRead_from_array (apers, AperData);

        if (!strcmp (AperData[0].name, refpixel))
          {
            if (AperData[0].data != NULL)
              {
                nbeam++;
              }
          }
        j = 0;
        while(AperData[j].name!=NULL){
          free(AperData[j++].data);
        }
      }
    return nbeam;
  }

/**
 * Function: nbeams_from_char_array
 * Read an Aperture File and returns an integer array
 * which lists the integer ID's of the beams found
 * in the given aperture.
 *
 * Parameters:
 * @param filename - name of Aperture File
 * @param num      - number of apertures found in Aperture File
 *
 * Returns:
 * @return beamIDs - the array with beam ID's
 */
gsl_vector_int * nbeams_from_char_array(char **apers, int num)
  {
    int i, j, nbeam = 0;
    char refpixel[MAXCHAR];

    gsl_vector_int *beamIDs;

    // allocate space for all possible beams
    beamIDs = gsl_vector_int_alloc(MAX_BEAMS);

    // set all integer ID's to the default value -1
    gsl_vector_int_set_all(beamIDs, -1);

    for (i = 0; i < MAX_BEAMS; i++)
      {
        struct CfgStrings AperData[] =
          {
            { refpixel, NULL },
            { NULL, NULL } /* array terminator. REQUIRED !!! */
          };
        sprintf(refpixel, "REFPIXEL%d%c", num, BEAM (i));
        CfgRead_from_array(apers, AperData);

        if (!strcmp(AperData[0].name, refpixel))
          {
            if (AperData[0].data != NULL)
              {
                // set the integer array to the ID
                gsl_vector_int_set(beamIDs, nbeam, i);
                nbeam++;
              }
          }

        j = 0;
        while (AperData[j].name!=NULL)
          free(AperData[j++].data);
      }
    return beamIDs;
  }


/**
 * Function: get_beam_from_char_array
 * Read a single beam from a single aperture from a char array
 * and set up the content of an object structure.
 *
 * Parameters:
 * @param filename - the name of the aperture file.
 * @param aperID   - the int ID of the aperture.
 * @param beamID   - the int ID of the beam.
 * @param b        - a beam object structure to contain the beam description
 *
 * Returns:
 * @return 0/1     - 0 if beam was not found, 1 if found
 */
int get_beam_from_char_array(char **aper, int aperID, int beamID, beam * b)
  {
    char refpixel[MAXCHAR], corners[MAXCHAR], curve[MAXCHAR];
    char width[MAXCHAR], orient[MAXCHAR];
    char slitgeom[MAXCHAR];
    char beamname[MAXCHAR], awidth[MAXCHAR], bwidth[MAXCHAR];
    char aorient[MAXCHAR], flux[MAXCHAR], modspec[MAXCHAR];
    char modimage[MAXCHAR], ignore[MAXCHAR];
    int ix;
    gsl_vector *v;
    gsl_vector *fvals;
    //beam *b = malloc(sizeof(beam));

    struct CfgStrings AperData[] =
      {
        { NULL, NULL },
        { NULL, NULL },
        { NULL, NULL },
        { NULL, NULL },
        { NULL, NULL },
        { NULL, NULL },
        { NULL, NULL },
        { NULL, NULL },
        { NULL, NULL },
        { NULL, NULL },
        { NULL, NULL },
        { NULL, NULL },
        { NULL, NULL },
        { NULL, NULL } /* array terminator. REQUIRED !!! */
      };

    AperData[0].name  = refpixel;
    AperData[1].name  = corners;
    AperData[2].name  = curve;
    AperData[3].name  = width;
    AperData[4].name  = orient;
    AperData[5].name  = slitgeom;
    AperData[6].name  = modspec;
    AperData[7].name  = modimage;
    AperData[8].name  = awidth;
    AperData[9].name  = bwidth;
    AperData[10].name = aorient;
    AperData[11].name = flux;
    AperData[12].name = ignore;

    sprintf(beamname, "%d%c", aperID, BEAM (beamID));

    sprintf(refpixel, "REFPIXEL%s", beamname);
    sprintf(corners,  "CORNERS%s",  beamname);
    sprintf(curve,    "CURVE%s",    beamname);
    sprintf(width,    "WIDTH%s",    beamname);
    sprintf(orient,   "ORIENT%s",   beamname);
    sprintf(slitgeom, "SLITGEOM%s", beamname);
    sprintf(modspec,  "MODSPEC%s",  beamname);
    sprintf(modimage, "MODIMAGE%s", beamname);
    sprintf(awidth,   "AWIDTH%s",   beamname);
    sprintf(bwidth,   "BWIDTH%s",   beamname);
    sprintf(aorient,  "AORIENT%s",  beamname);
    sprintf(flux,     "FLUX%s",     beamname);
    sprintf(ignore,   "IGNORE%s",   beamname);

    // read in the structure
    CfgRead_from_array(aper, AperData);

    for (ix = 0; ix < 13; ix++)
      {
        // get the ignore flag
        if (!strcmp(AperData[ix].name, ignore))
          b->ignore = atoi(AperData[ix].data);

        /* Get the reference point */
        if (!strcmp(AperData[ix].name, refpixel))
          {
            if (AperData[ix].data == NULL)
              {
                /* Whether the REFPIXEL entry is found or not determine whether this beam exists or not */
                fprintf(stderr, "%s exit!\n", beamname);
                return 0;
              }
            b->ID = beamID;
            v = string_to_gsl_array(AperData[ix].data);
            if (v->size != 2)
              {
                aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                "%s definition should contain 2 points.\n"
                "\"%s\" only appears to contain %i.",
                refpixel, AperData[ix].data, v->size);
              }
            b->refpoint.x = gsl_vector_get(v, 0);
            b->refpoint.y = gsl_vector_get(v, 1);
            gsl_vector_free(v);
          }

        /* Get the bounding box */
        if (!strcmp(AperData[ix].name, corners))
          {
            if (AperData[ix].data == NULL)
              aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
              "%s not found\n", corners);
            v = string_to_gsl_array(AperData[ix].data);
            if (v->size != 8)
              {
                aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                "%s definition must contain 8 (4x2) points.\n"
                "\"%s\" only appears to contain %i.",
                corners, AperData[ix].data, v->size);
              }
            b->corners[0].x = gsl_vector_get(v, 0);
            b->corners[0].y = gsl_vector_get(v, 1);
            b->corners[1].x = gsl_vector_get(v, 2);
            b->corners[1].y = gsl_vector_get(v, 3);
            b->corners[2].x = gsl_vector_get(v, 4);
            b->corners[2].y = gsl_vector_get(v, 5);
            b->corners[3].x = gsl_vector_get(v, 6);
            b->corners[3].y = gsl_vector_get(v, 7);
            gsl_vector_free(v);
          }

        /* Get and set up the trace */
        if (!strcmp(AperData[ix].name, curve))
          {
            if (AperData[ix].data == NULL)
              aXe_message(aXe_M_FATAL, __FILE__, __LINE__, "%s not found\n", curve);

            // convert the string to an array
            v = string_to_gsl_array(AperData[ix].data);
              {
                int ii;
                gsl_vector *vv;

                vv = gsl_vector_alloc(v->size -1);
                for (ii=0; ii<(int)v->size -1; ii++)
                  {
                    gsl_vector_set(vv, ii, gsl_vector_get(v, ii+1));
                  }
                b->spec_trace = create_polyN(vv);
                gsl_vector_free(vv);
              }
            gsl_vector_free(v);
          }

        /* Get the object orientation */
        if (!strcmp(AperData[ix].name, orient))
          {
            if (AperData[ix].data == NULL)
              aXe_message(aXe_M_FATAL, __FILE__, __LINE__, "%s not found in %s\n", orient);
            v = string_to_gsl_array(AperData[ix].data);
            if (v->size != 1)
              {
                aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                "%s definition must contain 1 angle (degrees).\n"
                "\"%s\" appears to contain %i.",
                orient, AperData[ix].data, v->size);
              }
            /* Convert from SeXtractor angle reference frame to aXe's */
            b->orient = (180 + gsl_vector_get(v, 0)) / 180. * M_PI;
            while (b->orient > M_PI)
              b->orient = b->orient - M_PI;
            gsl_vector_free(v);
          }

        /* Get the object width */
        if (!strcmp(AperData[ix].name, width))
          {
            if (AperData[ix].data == NULL)
              aXe_message(aXe_M_FATAL, __FILE__, __LINE__, "%s not found\n", width);
            v = string_to_gsl_array(AperData[ix].data);
            if (v->size != 1)
              {
                aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                "%s definition must contain 1 width (pixel).\n"
                "\"%s\" appears to contain %i.", width,
                AperData[ix].data, v->size);
              }
            b->width = gsl_vector_get(v, 0);
            gsl_vector_free(v);
          }

        // get the slit geometry
        if (!strcmp(AperData[ix].name, slitgeom))
          {
            // complain if no data is found
            //if (AperData[ix].data == NULL)
            //  aXe_message(aXe_M_FATAL, __FILE__, __LINE__, "%s not found\n", slitgeom);

            if (AperData[ix].data != NULL)
              {
                // convert the data string to a number vector
                v = string_to_gsl_array(AperData[ix].data);

                // check that you get four numbers
                if (v->size != 4)
                  aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                      "%s definition must contain 1 width (pixel).\n"
                      "\"%s\" appears to contain %i.", width,
                      AperData[ix].data, v->size);

                // transfer the data
                b->slitgeom[0] = gsl_vector_get(v, 0);
                b->slitgeom[1] = (gsl_vector_get(v, 1) + 180.0) / 180.0 * M_PI;
                b->slitgeom[2] = gsl_vector_get(v, 2);
                b->slitgeom[3] = gsl_vector_get(v, 3);

                // free memory
                gsl_vector_free(v);
              }
          }

        // get the spectral model index
        if (!strcmp(AperData[ix].name, modspec))
          {
            if (AperData[ix].data != NULL)
              b->modspec = atoi(AperData[ix].data);
          }

        // get the object shape index
        if (!strcmp(AperData[ix].name, modimage))
          {
            if (AperData[ix].data != NULL)
              b->modimage = atoi(AperData[ix].data);
          }

        /* Get the maximum object width */
        if (!strcmp(AperData[ix].name, awidth))
          {
            if (AperData[ix].data != NULL)
              {
                //      fprintf(stderr, "Maximum object is %s.\n", AperData[ix].data);
                v = string_to_gsl_array(AperData[ix].data);
                if (v->size != 1)
                  {
                    aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                    "%s definition must contain 1 width (pixel).\n"
                    "\"%s\" appears to contain %i.", width,
                    AperData[ix].data, v->size);
                  }
                b->awidth = gsl_vector_get(v, 0);
                gsl_vector_free(v);
              }
            else
              {
                // set the width to a default
                // value
                b->awidth = -1.0;
              }
          }

        /* Get the minimum object width */
        if (!strcmp(AperData[ix].name, bwidth))
          {
            if (AperData[ix].data != NULL)
              {
                //      fprintf(stderr, "Minimum object is %s.\n", AperData[ix].data);
                v = string_to_gsl_array(AperData[ix].data);
                if (v->size != 1)
                  {
                    aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                    "%s definition must contain 1 width (pixel).\n"
                    "\"%s\" appears to contain %i.", width,
                    AperData[ix].data, v->size);
                  }
                b->bwidth = gsl_vector_get(v, 0);
                gsl_vector_free(v);
              }
            else
              {
                // set the width to a default
                // value
                b->bwidth = -1.0;
              }
          }

        /* Get the orientation for maximum object width */
        if (!strcmp(AperData[ix].name, aorient))
          {
            if (AperData[ix].data != NULL)
              {
                //      fprintf(stderr, "Orienation for awidth is %s.\n", AperData[ix].data);
                v = string_to_gsl_array(AperData[ix].data);
                if (v->size != 1)
                  {
                    aXe_message(aXe_M_FATAL, __FILE__, __LINE__,
                    "%s definition must contain 1 width (pixel).\n"
                    "\"%s\" appears to contain %i.", width,
                    AperData[ix].data, v->size);
                  }
                //  Convert from SeXtractor angle reference frame to aXe's
                b->aorient = (180 + gsl_vector_get(v, 0)) / 180. * M_PI;
                while (b->aorient > M_PI)
                  b->aorient = b->aorient - M_PI;
                gsl_vector_free(v);
              }
          }

        /* Get the flux values */
        if (!strcmp(AperData[ix].name, flux))
          {
            if (AperData[ix].data != NULL)
              {
                //      fprintf(stderr, "Flux values are: %s.\n", AperData[ix].data);
                fvals = string_to_gsl_array(AperData[ix].data);
                b->flux = fvals;
                //      gsl_vector_free (v);
              }
          }

        /* Set up the bounding box */
        quad_to_bbox(b->corners, b->bbox, b->bbox + 1);

      }
    ix = 0;
    while (AperData[ix].name!=NULL)
      {
        free(AperData[ix++].data);
      }
    return 1;
  }


/**
 * Function: get_aperture_from_char_array
 * Read a specific aperture number from an array of strings
 *
 * Parameters:
 * @param filename - the name of the aperture file.
 * @param aperID   - the number of the aperture.
 *
 * Reaturns:
 * @return ob - an object type containing the aperture description and
 *              all the  included beams. Returns NULL if aperture was
 *              not found
 */
object *
get_aperture_from_char_array (char **aper, int aperID)
{
  int nbeams, ibeam = 0, beamID;
  object *ob = malloc (sizeof (object));
  gsl_vector_int *beamlist;

  ob->ID = aperID;
  //  nbeams = nbeams_from_char_array2 (aper, aperID);
  // get the array with the integer ID's for the beams
  beamlist = nbeams_from_char_array(aper, aperID);

  // derive the number of beams by checking
  // the integer array for non default values
  nbeams=0;
  while( nbeams < (int)beamlist->size && gsl_vector_int_get (beamlist, nbeams) > -1)
    nbeams++;

  // set some default values for struct
  // items which may NOT be overwritten with
  // OAF content
  for (beamID = 0; beamID < nbeams; beamID++)
    {
      ob->beams[beamID].flux        = NULL;
      ob->beams[beamID].slitgeom[0] = -1.0;
      ob->beams[beamID].slitgeom[1] = -1.0;
      ob->beams[beamID].slitgeom[2] = -1.0;
      ob->beams[beamID].slitgeom[3] = -1.0;
      ob->beams[beamID].modspec     = -1;
      ob->beams[beamID].modimage    = -1;
    }

  // go over all beams
  for (beamID = 0; beamID < nbeams; beamID++)
    // convert the current beam to a beam structure
    if (get_beam_from_char_array(aper, aperID, gsl_vector_int_get(beamlist,beamID),
                                   &(ob->beams[ibeam])))
      // enhance the counter
      ibeam++;

  // set the beam number
  ob->nbeams = ibeam;

  // free the integer ID array
  gsl_vector_int_free(beamlist);

  // return the object structure
  return ob;
}


/**
 * Function: file_to_object_list_seq
 * Read all apertures from n aperture file and return an array of object
 * pointers. This function use a sequential read of the aperture file which
 * requires apertures to be clearly separated by APERTURE #/APERTURE END
 * keys in the Aperture File. It is faster than the more generic
 * file_to_object_list function. This function relies on the
 * from_char-array() functions.
 *
 * Parameters:
 * @param filename - the name of the aperture file.
 * @param obs      - a pointer to an existing non NULL observation structure
 *
 * Returns:
 * @return oblist - a newly allocated array of object pointers
 */
object ** file_to_object_list_seq(char filename[], observation * obs)
  {
    gsl_vector_int *l;
    object **oblist;
    int i, j, aperID;
    char **aper;
    FILE *input;
    l = aper_file_aperlist(filename);
    if (l==NULL)
      return NULL; // no aperture was found
    /* Allocate memory for enough objects */
    oblist = (object **) malloc((l->size + 1) * sizeof(object *));

    // open the file;
    // report any problems
    input=fopen(filename, "r");
    if (input==NULL)
      aXe_message(aXe_M_FATAL, __FILE__, __LINE__, "Could not open %s", filename);

    for (i = 0; i < (int)l->size; i++)
      {
        aperID = gsl_vector_int_get(l, i);
        aper = return_next_aperture(input);
        oblist[i] = get_aperture_from_char_array(aper, aperID);
        oblist[i]->grism_obs = obs;

        // clean up memory allocated
        for (j=0; j<MAX_BEAMS*APER_MAXLINE; j++)
          free(aper[j]);
        free(aper);
      }

    oblist[i++] = NULL;
    fclose(input);
    gsl_vector_int_free(l);

    // return the object list
    return oblist;
  }

/**
 * Function: return_next_aperture
 * this function allocates and returns an array of strings containing
 * the content of the next APERTURE in an already opened aperture file.
 *
 * Parameters:
 * @param input - a pointer to an opened Aperture File
 *
 * Returns:
 * @return aper - a pointer to a newly allocated array of strings
 */
char** return_next_aperture(FILE *input)
  {
    char *WorkPtr;
    char Buffer[BUFFERSIZE], Buffer2[BUFFERSIZE];
    char *key;
    int new_aper=0;
    char **aper=NULL;
    int nlines= MAX_BEAMS*APER_MAXLINE, i=0;

    /* Allocate enough room to allow for MAXBEAM in returned string array */

    aper = (char **) malloc(sizeof(char *)*nlines);
    for (i=0; i<nlines; i++)
      aper[i] = malloc(sizeof(char)*BUFFERSIZE);

    i=0;
    while (NULL != fgets (Buffer, BUFFERSIZE, input))
      {
        strcpy(Buffer2,Buffer);
        /* clip off optional comment tail indicated by a semi-colon */
        if (NULL != (WorkPtr = strchr (Buffer, ';')))
          *WorkPtr = '\0';
        else
          WorkPtr = Buffer + strlen (Buffer);

        /* clip off trailing and leading white space*/
        WorkPtr--;
        while (isspace ((int) *WorkPtr) && WorkPtr >= Buffer)
          *WorkPtr-- = '\0';

        WorkPtr = Buffer;
        while (isspace ((int) *WorkPtr))
          WorkPtr++;

        if (0 == strlen (WorkPtr))
          continue;
        key = strtok (WorkPtr, " =");

        if (!strcmp(key,"APERTURE"))
          new_aper++;

        sprintf(aper[i++], "%s", Buffer2);

        if (new_aper==2)
          {
            sprintf(aper[i++],"NULL");
            return aper;
          }
        }
    return aper;
  }


/**
 * Function: find_object_in_object_list
 * This function returns the array index of the object in a NULL terminated
 * array of objects which corresponds to the passed ID
 *
 * Parameters:
 * @param oblist - a NULL terminated array of objects
 * @param ID     - the numeric ID of the object/aperture to look for
 *
 * Returns:
 * @return i/-1  - the array index to the object/aperture with the wanted ID
 */
int
find_object_in_object_list(object **oblist, const int ID)
{
  int i;

  i=0;
  while (oblist != NULL && oblist[i] != NULL)
    {
      if (oblist[i]->ID == ID) return i;
      i++;
    }

  //    return NULL;
  // this is critical and not throughout testet!
  return -1;
}


/**
 * Function: find_beam_in_object_list
 * This function returns the array index of the object in a NULL terminated
 * array of objects which corresponds to the passed ID
 *
 * Parameters:
 * @param oblist - a NULL terminated array of objects
 * @param objID  - the numeric ID of the object/aperture to look for
 * @param beamID - the numeric ID of the beam to look for
 *
 * Returns:
 * @return actbeam  - the beam identitified
 */
beam
find_beam_in_object_list(object **oblist, const int objID,
                         const int beamID)
{
  beam actbeam;

  int obj_index=0;
  int j=0;

  // set the beam ID to -1 to identify
  // failed identification
  actbeam.ID = -1;

  // identify the correct index in the object list
  obj_index = find_object_in_object_list(oblist, objID);

  //
  if (obj_index > -1)
    {
      // go over all beams in the matchin object
      for (j=0; j < oblist[obj_index]->nbeams; j++)
        {

          // search for a matching beam ID
          if (oblist[obj_index]->beams[j].ID == beamID)
            // copy the matching beam
            actbeam = oblist[obj_index]->beams[j];
        }
    }

  // return the beam
  return actbeam;
}

/**
 * Function: find_beamptr_in_object_list
 * This function returns the array index of the object in a NULL terminated
 * array of objects which corresponds to the passed ID
 *
 * Parameters:
 * @param oblist - a NULL terminated array of objects
 * @param objID  - the numeric ID of the object/aperture to look for
 * @param beamID - the numeric ID of the beam to look for
 *
 * Returns:
 * @return actbeam  - the beam identitified
 */
beam *
find_beamptr_in_object_list(object **oblist, const int objID,
                            const int beamID)
{
  beam *actbeam=NULL;

  int obj_index=0;
  int j=0;

  // set the beam ID to -1 to identify
  // failed identification
  //actbeam->ID = -1;

  // identify the correct index in the object list
  obj_index = find_object_in_object_list(oblist, objID);
  //fprintf(stderr, "found index: %i\n", obj_index);
  //
  if (obj_index > -1)
    {
      // go over all beams in the matchin object
      for (j=0; j < oblist[obj_index]->nbeams; j++)
        {

          // search for a matching beam ID
          if (oblist[obj_index]->beams[j].ID == beamID)
            // copy the matching beam
            actbeam = &oblist[obj_index]->beams[j];
        }
    }

  // return the beam
  return actbeam;
}

void
refurbish_object_list(object **oblist, const int new_default,
                       const int old_value, const int new_value)
{
  //object *act_obj;
  beam   *act_beam;

  int j=0;
  int i=0;

  i=0;
  while (oblist != NULL && oblist[i] != NULL)
    {
      // go over all beams in the matchin object
      for (j=0; j < oblist[i]->nbeams; j++)
        {
          act_beam = &oblist[i]->beams[j];
          if (act_beam->ignore == old_value)
            act_beam->ignore = new_value;
          else
            act_beam->ignore = new_default;
        }
      i++;
    }
}

/**
 * Function: object_list_size
 * This function returns the number of element in an array
 * of objects.
 *
 * Parameters:
 * @param oblist - a NULL terminated array of objects
 *
 * Returns:
 * @return i     - number of elements
 */
int object_list_size(object **oblist)
{
    int i=0;

    if (oblist==NULL) return 0;

    while (oblist[i] != NULL)
      i++;

    return i;
}

/**
 * Function: get_beamspec_size
 * The function computes the numbers of beams that
 * are supposed to modelled. This is done by checking
 * the ignore flag of each beam in an object list.
 * Some of the beams might be completely outside of
 * the image are. So the number delivered here is only
 * an upper limit.
 *
 * Parameter:
 * @param  oblist   - the list of objects
 *
 * Returns:
 * @return specsize - the number of objects to be modeled
 */
int get_beamspec_size(object **oblist)
{
  int specsize=0;
  int objectsize=0;
  int i=0;
  int j=0;

  // determine the number of objects
  objectsize = object_list_size(oblist);

  // go over each object
  for (i=0; i< objectsize; i++)
    {
      // go over each beam
      for (j=0; j < oblist[i]->nbeams; j++)
        {
          // check the ignore flag
          // sum up the counter if permitted
          if (oblist[i]->beams[j].ignore != 1)
            specsize++;
        }
    }

  // return the counter
  return specsize;
}
