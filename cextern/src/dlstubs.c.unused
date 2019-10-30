/**
 * File: dlstubs.c
 * Some empty functions needed at some time to statically
 * link in Solaris.
 * @author  Martin Kuemmel
 * @package dlstubs
 * @version $Revision: 1.2 $
 * @date    $Date: 2010-06-15 09:48:34 $
 */
#include <sys/types.h>
#include <dlfcn.h>

/* dl*() stub routines for static compilation.  Prepared from
   /usr/include/dlfcn.h by Hal Pomeranz <hal@deer-run.com> */

void *dlopen(const char *str, int x) {}
void *dlsym(void *ptr, const char *str) {}
int dlclose(void *ptr) {}
char *dlerror() {}
void *dlmopen(Lmid_t a, const char *str, int x) {}
int dladdr(void *ptr1, Dl_info *ptr2) {}
int dldump(const char *str1, const char *str2, int x) {}
int dlinfo(void *ptr1, int x, void *ptr2) {}

void *_dlopen(const char *str, int x) {}
void *_dlsym(void *ptr, const char *str) {}
int _dlclose(void *ptr) {}
char *_dlerror() {}
void *_dlmopen(Lmid_t a, const char *str, int x) {}
int _dladdr(void *ptr1, Dl_info *ptr2) {}
int _dldump(const char *str1, const char *str2, int x) {}
int _dlinfo(void *ptr1, int x, void *ptr2) {}
