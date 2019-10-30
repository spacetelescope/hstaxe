#include "inima_utils.h"
#include "spc_cfg.h"



/**
    Utility function to parse a filename and populate a CfgStrings structure
    
    @param Filename a pointer to a string containing the name of a configuration file
    @param CfgInfo a pointer to an existing CfgStrings structure
    @return 

    @see CfgStrings

*/
int
CfgRead (char *Filename, struct CfgStrings *CfgInfo)
{
  static char Buffer[BUFFERSIZE];
  char *WorkPtr;
  char *CfgName;
  char *CfgData;
  struct CfgStrings *Cfg;
  FILE *CfgFile;
  
  CfgFile = fopen (Filename, "r");
  if (NULL == CfgFile)
    {
      aXe_message (aXe_M_FATAL, __FILE__, __LINE__, "Could not open %s!\n",
		   Filename); 
    }

  //  fprintf(stderr, "%s: %s\n", Filename, CfgInfo->name);

  while (NULL != fgets (Buffer, BUFFERSIZE, CfgFile))
    {
      // check whether it is a valid
      // line, or whether it can be discarded
      if (!is_valid_inima_line(Buffer))
      	continue;
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
      WorkPtr = NULL;

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
	  //strlwr( CfgName );
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

	  /* look for matching 'name'
	   */
	  Cfg = CfgInfo;
	  while (NULL != Cfg->name
		 && 0 != strcmp (Cfg->name, CfgName))
	    Cfg++;

	  /* duplicate the data if the name is found.
	   */
	  //      fprintf(stderr, "%s: %s: %s %i\n", Filename, CfgInfo->name, Buffer, strlen(Buffer));
	  //      fprintf(stderr, "%s: %s:\n", Cfg->name, CfgData);
	  if (NULL != Cfg->name)
	    {
	      if (NULL == CfgData)
		{
		  fclose (CfgFile);
		  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			       "Error. No data read for %s. Check format of configuration file.\n",
			       Cfg->name);
		}
	      Cfg->data = strdup (CfgData);	/* strdup is not ANSI    */
	      /* memory leaks if Cfg->data */
	      /* is malloc'ed already      */
	      
	    }		/* undetected error on failure should not be a problem  */
		   /* as configuration reading should be done early.       */
	}			/* but test and handle it anyway ...                    */
    }
  fclose (CfgFile);
  return NO_PROBLEMS;
}


/**
    Utility function to parse a NULL terminated array of strings and populate a CfgStrings structure
    
    @param Filename a pointer to a string containing the name of a confiuration file
    @param CfgInfo a pointer to an existing CfgStrings structure
    @return 

    @see CfgStrings

*/
int
CfgRead_from_array (char **arr, struct CfgStrings *CfgInfo)
{
  char Buffer[BUFFERSIZE];
  char *WorkPtr;
  char *CfgName;
  char *CfgData;
  struct CfgStrings *Cfg;
  int i=0;
  
  while (NULL != arr[i] && strcmp(arr[i],"NULL\0") != 0)
    {
      sprintf(Buffer,"%s",arr[i]);
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
	  //strlwr( CfgName );
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
	  
	  /* look for matching 'name'
	   */
	  Cfg = CfgInfo;
	  while (NULL != Cfg->name
		 && 0 != strcmp (Cfg->name, CfgName))
	    Cfg++;
	  
	  /* duplicate the data if the name is found.
	   */
	  if (NULL != Cfg->name)
	    {
	      Cfg->data = strdup (CfgData);	/* strdup is not ANSI    */
	      //	      strcpy(Cfg->data,CfgData);	/* strdup is not ANSI    */
	      /* memory leaks if Cfg->data */
	      /* is malloc'ed already      */
	      if (NULL == Cfg->data)
		{
		  aXe_message (aXe_M_FATAL, __FILE__, __LINE__,
			       "Memory error. No data read");
		}
	    }		/* undetected error on failure should not be a problem  */
	  /* as configuration reading should be done early.       */
	}			/* but test and handle it anyway ...                    */
      i++;
    }
  return NO_PROBLEMS;
}
