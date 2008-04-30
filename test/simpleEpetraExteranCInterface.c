#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "FEpetra_Map.h"
#include "FEpetra_Vector.h"


int main( int argc, char* argv[] )
{
  
  /*
   * Data declarations (how old-school is this!)
   */

  int numGlobalElements, numGlobalElements_rtn;
  MapID mapID;
  VectorID xID, bID;
  double bnorm, xnorm, expected_bnorm, expected_xnorm, bnorm_err, xnorm_err, err_tol;
  int success = 1;

  /*
   * Executable code
   */
  
  /* Create a map */
  numGlobalElements = 4;
  mapID = FEpetra_Map_Create(numGlobalElements);

  numGlobalElements_rtn = FEpetra_Map_NumGlobalElements(mapID);
  printf( "NumGlobalElements = %d\n", numGlobalElements_rtn );
  assert( numGlobalElements == numGlobalElements_rtn );
  
  /* Create vectors */
  xID = FEpetra_Vector_Create(mapID);
  bID = FEpetra_Vector_Create(mapID);

  /* Do some vector operations */
  FEpetra_Vector_PutScalar(bID, 2.0);
  FEpetra_Vector_Update(xID, 2.0, bID, 0.0); /* x = 2*b */

  bnorm = FEpetra_Vector_Norm2(bID);
  xnorm = FEpetra_Vector_Norm2(xID);

  printf( "2 norm of x = %f\n", xnorm );
  printf( "2 norm of b = %f\n", bnorm );

  /* Test the expected value */
  err_tol = 0.0;
  expected_bnorm = sqrt( 2.0 * 2.0 * numGlobalElements );
  expected_xnorm = sqrt( 4.0 * 4.0 * numGlobalElements );
  bnorm_err = fabs( expected_bnorm - bnorm ) / expected_bnorm;
  xnorm_err = fabs( expected_xnorm - xnorm ) / expected_xnorm;
  printf( "error in 2 norm of x = %f\n", bnorm_err );
  printf( "error in 2 norm of b = %f\n", xnorm_err );
  if (bnorm_err > err_tol) success = 0;
  if (xnorm_err > err_tol) success = 0;

  /* Clean up memory (in reverse order)! */
  FEpetra_Vector_Destroy(bID);
  FEpetra_Vector_Destroy(xID);
  FEpetra_Map_Destroy(mapID);

  /* This should throw an exception and print an error message! */
  /* FEpetra_Map_NumGlobalElements(mapID); */

  if (success)
    printf( "\nEnd Result: TEST PASSED\n" );
  else
    printf( "\nEnd Result: TEST FAILED\n" );
  
  return ( success ? 0 : 1 );

}
