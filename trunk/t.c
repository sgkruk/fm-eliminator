#include "setoper.h"
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

dd_boolean SetInputFile(FILE **f, dd_DataFileType fname)
{
  dd_boolean success=dd_FALSE;
  success=dd_FALSE;

  if ( ( *f = fopen(fname, "r") )!= NULL) {
    printf("input file %s is open\n", fname);
    success=dd_TRUE;
  }
  else{
    printf("The input file %s not found\n",fname);
  }
  return success;
}

dd_boolean SetWriteFile(FILE **f, dd_DataFileType fname)
{
  dd_boolean success=dd_FALSE;

  if ( (*f = fopen(fname, "w")) != NULL){
    printf("output file %s is open\n",fname);
    success=dd_TRUE;
  }
  else{
    printf("The output file %s cannot be opened\n",fname);
  }
  return success;
}


int main(int argc, char *argv[])
{
  dd_MatrixPtr M=NULL,M1=NULL;
  dd_colrange j,s,t,d;
  dd_ErrorType err=dd_NoError;
  dd_rowset redset,impl_linset;
  dd_colset delset;
  dd_rowindex newpos;
  mytype val;
  dd_DataFileType inputfile;
  FILE *reading=NULL;

  dd_set_global_constants();  /* First, this must be called. */

  dd_init(val);
  




  dd_FreeMatrix(M1);
  set_free(delset);
  set_free(redset);
  set_free(impl_linset);
  free(newpos);

_L99:;
  /* if (err!=dd_NoError) dd_WriteErrorMessages(stderr,err); */
  dd_free_global_constants();  /* At the end, this should be called. */
  return 0;
}


/* end of projection.c */

