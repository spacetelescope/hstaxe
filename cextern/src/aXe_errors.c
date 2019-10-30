/* 
 */

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include "aXe_errors.h"

/**
  Report an error message or warning and exit if error is fatal.

  @param level severety of the message, one of <ul><li>aXe_M_FATAL (fatal error,
    exit),</li>
    <li>aXe_M_ERROR (an error that will invalidate the result but doesn't
    preclude further processing within the program), </li>
    <li>aXe_M_WARN1 (a problem that will most likely invalidate the result),
    </li>
    <li>aXe_M_WARN2 (a problem that may invalidate the result),</li>
    <li>aXe_M_WARN3 (a problem that might be looked into, but probably all is 
    ok),</li>
    <li>aXe_M_WARN4 (informational messages only).</li></ul>
  @param w_file the module the error occurred in.
  @param w_line the line the error occurred in.
  @param message the error message -- this works like printf
*/
void
aXe_message (const int level, const char *const w_file, const int w_line,
	     const char *const msg, ...)
/* This is based on fvwm's fvwm_msg */
{
     char *typestr;
     va_list args;

     switch (level)
       {
       case aXe_M_FATAL:
	    typestr = "Fatal:";
	    break;
       case aXe_M_ERROR:
	    typestr = "Error:";
	    break;
       case aXe_M_WARN1:
       case aXe_M_WARN2:
       case aXe_M_WARN3:
       case aXe_M_WARN4:
	    typestr = "Warning:";
	    break;
       default:
	    typestr = "Oops:";
	    break;
       }

     va_start (args, msg);

     fprintf (stderr, "aXe (%s, %d): %s ", w_file, w_line, typestr);
     vfprintf (stderr, msg, args);
     //fprintf (stderr, "\n");

     va_end (args);

     if (level == aXe_M_FATAL)
	  exit (-1);
}
