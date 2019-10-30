#ifndef _SPC_CFG_H
#define _SPC_CFG_H

#define LINE_LEN_MAX    1280	/* actual max line length  */
#define BUFFERSIZE      LINE_LEN_MAX +2	/* ... including \n and \0 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "aXe_errors.h"
#include "spc_utils.h"


enum RetVal
{
     NO_PROBLEMS,
     ERR_FOPEN, 
     ERR_MEM,
};

struct CfgStrings
{
     char *name;
     char *data;
};

extern int
CfgRead (char *Filename, struct CfgStrings *CfgInfo);

extern int
CfgRead_from_array (char **arr, struct CfgStrings *CfgInfo);
#endif
