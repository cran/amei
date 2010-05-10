#include "matrix.h"
#include "ssdiff.h"
#include <math.h>

/********************************************************/
/* ssdiff  */
/* PURPOSE: finds the sum of the square difference  */
/*				of 2 matrices.  */
/* INPUT: A and B and their dims M1, M2  */
/* OUTPUT: sum  */
/********************************************************/

double ssdiff(double **A, double **B, unsigned int M1,
	      unsigned int M2)
{
  double sum;	
  double **C;
  unsigned int  i, j;
  
  C = matrixdiff(A,B,M1,M2);
  
  sum = 0;
  
  for(i=0; i<M1; i++){
    for(j=0; j<M2; j++){
      sum= sum + C[i][j] * C[i][j];
    }
  }
  
  delete_matrix(C);
  
  return(sqrt(sum)/(M1*M2));
}

/********************************************************/
/* matrixdiff  */
/* PURPOSE: takes the difference of 2 matrices  */
/* INPUT: A, B, M1 rows, and M2 columns  */
/* OUTPUT: out (matrix of differences)  */
/********************************************************/

double** matrixdiff(double **A, double **B, unsigned int M1, 
		    unsigned int M2)
{
  double** out;
  unsigned int  i, j;
  
  out = new_zero_matrix(M1,M2);
  
  for(i=0; i<M1; i++){
    for(j=0; j<M2; j++){
      out[i][j] = A[i][j] - B[i][j];
    }
  }
  
  return(out);
  
}
