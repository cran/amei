#include "matrix.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TRUE 1
#define FALSE 0

int ivequal(int* v1, int* v2, int n)
{
  int k;
  for(k=0;k<n;k++)
    if(v1[k]!=v2[k])
      return 0;
  return 1;
}

int dvequal(double* v1, double* v2, int n)
{
  int k;
  for(k=0;k<n;k++)
    if(v1[k]!=v2[k])
      return 0;
  return 1;
}

double** diag(double d, int n)
{
  int k;
  double **m;
  m = new_2ddarray(n,n);
  zerov(*m,n*n);
  for(k=0;k<n;k++)
    m[k][k] = d;
  return m;
}

void extractSubMatrix(double** A,int dim,int* ind,double** Asub)
{
  int i,j;
  for(i=0;i<dim;i++)
    for(j=0;j<dim;j++)
      Asub[i][j] = A[ind[i]][ind[j]];
 }

void insertSubMatrix(double** Asub, int dim, int* ind, double** A)
{
  int i,j;
  for(i=0;i<dim;i++)
    for(j=0;j<dim;j++)
      A[ind[i]][ind[j]] = Asub[i][j];
}

float** new_2dfarray(int d1, int d2)
{
  float **a;
  int i;
  a = (float**) malloc(sizeof(float*)*d1);
  a[0] = (float*) malloc(sizeof(float)*d1*d2);
  for(i=1;i<d1;i++) a[i] = a[i-1]+d2;
  return a;
}

void rm_2dfarray(float** x)
{
  free(x[0]);
  free(x);
  return;
}

int** new_2diarray(int d1, int d2)
{
  int **a;
  int i;
  a = (int**) malloc(sizeof(int*)*d1);
  a[0] = (int*) malloc(sizeof(int)*d1*d2);
  for(i=1;i<d1;i++) a[i] = a[i-1]+d2;
  return a;
}

void rm_2diarray(int** x)
{
  free(x[0]);
  free(x);
  return;
}

double** new_2ddarray(int d1, int d2)
{
  double **a;
  int i;
  a = (double**) malloc(sizeof(double*)*d1);
  a[0] = (double*) malloc(sizeof(double)*d1*d2);
  for(i=1;i<d1;i++) a[i] = a[i-1]+d2;
  return a;
}

void rm_2ddarray(double** x)
{
  free(x[0]);
  free(x);
  return;
}

float** new_2dfptrs(int d1, int d2, float* x)
{
  float **a;
  int i;
  a = (float**) malloc(sizeof(float*)*d1);
  a[0] = x;
  for(i=1;i<d1;i++) a[i] = a[i-1]+d2;
  return a;
}

void rm_2dfptrs(float** x)
{
  free(x);
  return;
}

int** new_2diptrs(int d1, int d2, int* x)
{
  int **a;
  int i;
  a = (int**) malloc(sizeof(int*)*d1);
  a[0] = x;
  for(i=1;i<d1;i++) a[i] = a[i-1]+d2;
  return a;
}

void rm_2diptrs(int** x)
{
  free(x);
  return;
}

double** new_2ddptrs(int d1, int d2, double* x)
{
  double **a;
  int i;
  a = (double**) malloc(sizeof(double*)*d1);
  a[0] = x;
  for(i=1;i<d1;i++) a[i] = a[i-1]+d2;
  return a;
}

void rm_2ddptrs(double** x)
{
  free(x);
  return;
}


/* 
   initialize a new array of floats with dimensions [d1] x [d2] x [d3]
*/
float*** new_3dfarray(int  d1, int d2, int d3)
{
  float*** a;
  int i,j;

  a = (float***) malloc(sizeof(float**)*d1);

  for(i=0;i<d1;i++)
    a[i] = (float**) malloc(sizeof(float*)*d2);

  a[0][0] = (float*) malloc(sizeof(float)*d1*d2*d3);
  for(i=1;i<d1;i++) a[i][0] = a[i-1][0] + (d2*d3);
  for(i=0;i<d1;i++) for(j=1;j<d2;j++) a[i][j] = a[i][j-1] + d3;
  
  return a;
}


/*
  deletes the array of floats initialized as above
*/
void rm_3dfarray(float*** myarray, int d1)
{
  int i;
  free(myarray[0][0]);
  for(i=0;i<d1;i++)
    free(myarray[i]);
  free(myarray);
  return;
}


/* 
   initialize a new array of ints with dimensions [d1] x [d2] x [d3]
*/
int*** new_3diarray(int  d1, int d2, int d3)
{
  int*** a;
  int i,j;

  a = (int***) malloc(sizeof(int**)*d1);

  for(i=0;i<d1;i++)
    a[i] = (int**) malloc(sizeof(int*)*d2);

  a[0][0] = (int*) malloc(sizeof(int)*d1*d2*d3);
  for(i=1;i<d1;i++) a[i][0] = a[i-1][0] + (d2*d3);
  for(i=0;i<d1;i++) for(j=1;j<d2;j++) a[i][j] = a[i][j-1] + d3;
  
  return a;
}


/*
  deletes the array of ints initialized as above
*/
void rm_3diarray(int*** myarray, int d1)
{
  int i;
  free(myarray[0][0]);
  for(i=0;i<d1;i++)
    free(myarray[i]);
  free(myarray);
  return;
}



/* 
   initialize a new array of doubles with dimensions [d1] x [d2] x [d3]
*/
double*** new_3ddarray(int  d1, int d2, int d3)
{
  double*** a;
  int i,j;

  a = (double***) malloc(sizeof(double**)*d1);

  for(i=0;i<d1;i++)
    a[i] = (double**) malloc(sizeof(double*)*d2);

  a[0][0] = (double*) malloc(sizeof(double)*d1*d2*d3);
  for(i=1;i<d1;i++) a[i][0] = a[i-1][0] + (d2*d3);
  for(i=0;i<d1;i++) for(j=1;j<d2;j++) a[i][j] = a[i][j-1] + d3;
  
  return a;
}

/*
  deletes the array of doubles initialized as above
*/
void rm_3ddarray(double*** myarray, int d1)
{
  int i;
  free(myarray[0][0]);
  for(i=0;i<d1;i++)
    free(myarray[i]);
  free(myarray);
  return;
}


/* 
   initializes the pointers to an array of floats with dimensions [d1] x [d2] x [d3]
   does NOT initialize memory for the actual floats, which are pointed to by x  
*/
float*** new_3dfptrs(int  d1, int d2, int d3, float* x)
{
  float*** a;
  int i,j;

  a = (float***) malloc(sizeof(float**)*d1);

  for(i=0;i<d1;i++)
    a[i] = (float**) malloc(sizeof(float*)*d2);

  a[0][0] = x;
  for(i=1;i<d1;i++) a[i][0] = a[i-1][0] + (d2*d3);
  for(i=0;i<d1;i++) for(j=1;j<d2;j++) a[i][j] = a[i][j-1] + d3;
  
  return a;
}

/*
  deletes the arrays of pointers initialized above
*/
void rm_3dfptrs(float*** myarray, int d1)
{
  int i;
  for(i=0;i<d1;i++)
    free(myarray[i]);
  free(myarray);
  return;
}

/* 
   initializes the pointers to an array of ints with dimensions [d1] x [d2] x [d3]
   does NOT initialize memory for the actual ints, which are pointed to by x 
*/
int*** new_3diptrs(int  d1, int d2, int d3, int* x)
{
  int*** a;
  int i,j;

  a = (int***) malloc(sizeof(int**)*d1);

  for(i=0;i<d1;i++)
    a[i] = (int**) malloc(sizeof(int*)*d2);

  a[0][0] = x;
  for(i=1;i<d1;i++) a[i][0] = a[i-1][0] + (d2*d3);
  for(i=0;i<d1;i++) for(j=1;j<d2;j++) a[i][j] = a[i][j-1] + d3;
  
  return a;
}

/*
  deletes the arrays of pointers initialized above
*/
void rm_3diptrs(int*** myarray, int d1)
{
  int i;
  for(i=0;i<d1;i++)
    free(myarray[i]);
  free(myarray);
  return;
}


/* 
   initializes the pointers to an array of doubles with dimensions [d1] x [d2] x [d3]
   does NOT initialize memory for the actual doubles, which are pointed to by x  
*/
double*** new_3ddptrs(int  d1, int d2, int d3, double* x)
{
  double*** a;
  int i,j;

  a = (double***) malloc(sizeof(double**)*d1);

  for(i=0;i<d1;i++)
    a[i] = (double**) malloc(sizeof(double*)*d2);

  a[0][0] = x;
  for(i=1;i<d1;i++) a[i][0] = a[i-1][0] + (d2*d3);
  for(i=0;i<d1;i++) for(j=1;j<d2;j++) a[i][j] = a[i][j-1] + d3;
  
  return a;
}

/*
  deletes the arrays of pointers initialized above
*/
void rm_3ddptrs(double*** myarray, int d1)
{
  int i;
  for(i=0;i<d1;i++)
    free(myarray[i]);
  free(myarray);
  return;
}


double**** new_4ddptrs(int d1, int d2, int d3, int d4, double* x)
{
  int n,i,j;
  double**** a;
  
  a = (double****) malloc(sizeof(double***)*d1);
  for(n=0;n<d1;n++)
    {
      a[n] = (double***) malloc(sizeof(double**)*d2);
      for(i=0;i<d2;i++)
	{
	  a[n][i] = (double**) malloc(sizeof(double*)*d3);
	}
    }
  a[0][0][0] = x;
  for(n=1;n<d1;n++) a[n][0][0] = a[n-1][0][0]+(d2*d3*d4);
  for(n=0;n<d1;n++) for(i=1;i<d2;i++) a[n][i][0] = a[n][i-1][0]+(d3*d4);
  for(n=0;n<d1;n++) for(i=0;i<d2;i++) for(j=1;j<d3;j++) a[n][i][j] = a[n][i][j-1] + d4;

  return a;
}

void rm_4ddptrs(double**** x, int d1, int d2)
{
  int n,i;
  for(n=0;n<d1;n++)
    {
      for(i=0;i<d2;i++)
	free(x[n][i]);
      free(x[n]);
    }
  free(x);
  return;
}

double**** new_4ddarray(int d1, int d2, int d3, int d4)
{
  int n,i,j;
  double**** a;
  
  a = (double****) malloc(sizeof(double***)*d1);
  for(n=0;n<d1;n++)
    {
      a[n] = (double***) malloc(sizeof(double**)*d2);
      for(i=0;i<d2;i++)
	{
	  a[n][i] = (double**) malloc(sizeof(double*)*d3);
	}
    }
  a[0][0][0] = (double*) malloc(sizeof(double)*d1*d2*d3*d4);
  for(n=1;n<d1;n++) a[n][0][0] = a[n-1][0][0]+(d2*d3*d4);
  for(n=0;n<d1;n++) for(i=1;i<d2;i++) a[n][i][0] = a[n][i-1][0]+(d3*d4);
  for(n=0;n<d1;n++) for(i=0;i<d2;i++) for(j=1;j<d3;j++) a[n][i][j] = a[n][i][j-1] + d4;

  return a;
}

void rm_4ddarray(double**** x, int d1, int d2)
{
  int n,i;
  
  free(x[0][0][0]);
  for(n=0;n<d1;n++)
    {
      for(i=0;i<d2;i++)
	free(x[n][i]);
      free(x[n]);
    }
  free(x);
  return;
}

void swap_3ddarray(double*** a, double*** b, int d1, int d2, int d3)
{
  double* tmp;
  int i,j;
  tmp = a[0][0];
  
  a[0][0] = b[0][0];
  b[0][0] = tmp;
  
  for(i=1;i<d1;i++)
    {
      a[i][0] = a[i-1][0]+(d2*d3);
      b[i][0] = b[i-1][0]+(d2*d3);
    }
  for(i=0;i<d1;i++)
    for(j=1;j<d2;j++)
      {
	a[i][j] = a[i][j-1]+d3;
	b[i][j] = b[i][j-1]+d3;
      }
  return;
}

/*
 * replace matrix with zeros
 */

void zero(double **M, unsigned int n1, unsigned int n2)
{
	unsigned int i, j;
	for(i=0; i<n1; i++) for(j=0; j<n2; j++) M[i][j] = 0;
}


/*
 * replace square matrix with identitiy 
 */

void id(double **M, unsigned int n)
{
	unsigned int i;
	zero(M, n, n);
	for(i=0; i<n; i++) M[i][i] = 1.0;
}


/*
 * same as new_matrix below, but for creating
 * n x n identity matrices
 */

double ** new_id_matrix(unsigned int n)
{
  unsigned int i;
  double** m = new_zero_matrix(n, n);
  for(i=0; i<n; i++) m[i][i] = 1.0;
  return m;
}


/*
 * same as new_matrix below, but zeros out the matrix
 */

double ** new_zero_matrix(unsigned int n1, unsigned int n2)
{
  unsigned int i, j;
  double **m = new_matrix(n1, n2);
  for(i=0; i<n1; i++) for(j=0; j<n2; j++) m[i][j] = 0.0;
  return m;
}


/*
 * create a new n1 x n2 matrix which is allocated like
 * and n1*n2 array, but can be alloced with [][]
 */
double ** new_matrix(unsigned int n1, unsigned int n2)
{
  int i;
  double **m;
  
  if(n1 == 0 || n2 == 0) return NULL;
  
  m = (double**) malloc(sizeof(double*) * n1);
  assert(m);
  m[0] = (double*) malloc(sizeof(double) * (n1*n2));
  assert(m[0]);
  
  for(i=1; i<n1; i++) m[i] = m[i-1] + n2;
  
  return m;
}


/*
 * create a new n1 x n2 matrix which is allocated like
 * and n1*n2 array, and copy the of n1 x n2 M into it.
 */

double ** new_dup_matrix(double** M, unsigned int n1, unsigned int n2)
{
  double **m;
  
  if(n1 <= 0 || n2 <= 0) {
    assert(M == NULL);
    return NULL;
  }
  
  m = new_matrix(n1, n2);
  dup_matrix(m, M, n1, n2);
  return m;
}


/*
 * copy M2 to M1
 */

void dup_matrix(double** M1, double **M2, unsigned int n1, unsigned int n2)
{
  assert(M1 && M2);
  unsigned int i, j;
  for(i=0; i<n1; i++) for(j=0; j<n2; j++) M1[i][j] = M2[i][j];
}


/*
 * swap the pointers of M2 to M1, and vice-versa
 * (tries to aviod dup_matrix when unnecessary)
 */

void swap_matrix(double **M1, double **M2, unsigned int n1, unsigned int n2)
{
  unsigned int  i;
  double *temp;
  temp = M1[0];
  M1[0] = M2[0];
  M2[0] = temp;
  for(i=1; i<n1; i++) {
    M1[i] = M1[i-1] + n2;
    M2[i] = M2[i-1] + n2;
  }
}

/*
 * delete a matrix allocated as above
 */
void delete_matrix(double** m)
{
  assert(*m);
  free(*m);
  assert(m);
  free(m);
}

/*
 * print an n x col matrix allocated as above out an opened outfile.
 * actually, this routine can print any double** 
 */
void printMatrix(double **M, unsigned int n, unsigned int col, FILE *outfile)
{
  int i,j;
  
  for(i=0; i<n; i++) {
    for(j=0; j<col; j++) {
      if(j==col-1) fprintf(outfile, "%.6f\n", M[i][j]);
      else fprintf(outfile, "%.6f ", M[i][j]);
    }
  }
}

/*
 * allocate and return an array of length n with scale*1 at
 * each entry
 */
double* ones(unsigned int n, double scale)
{
  double *o;
  unsigned int i;
  o = (double*) malloc(sizeof(double) * n);
  assert(o);
  for(i=0; i<n; i++) o[i] = scale;
  return o;
}

/*
 * allocate and return an array containing
 * the seqence of doubles [from...to] with steps of
 * size by
 */
double* dseq(double from, double to, double by)
{
  unsigned int n,i;
  double *s = NULL;
  
  by = abs(by);
  
  if(from <= to) n = (unsigned int) (to - from)/abs(by) + 1;
  else n = (unsigned int) (from - to)/abs(by) + 1;
  
  if( n == 0 ) return NULL;
  
  s = (double*) malloc(sizeof(double) * n);
  assert(s);
  s[0] = from;
  for(i=1; i<n; i++) {
    s[i] = s[i-1] + by;
  }
  return s;
}

/*
 * allocate and return an array containing
 * the integer seqence [from...to]
 */

int* iseq(double from, double to)
{
  unsigned int n,i;
  int by;
  int *s = NULL;
  
  if(from <= to) {
    n = (unsigned int) (to - from) + 1;
    by = 1;
  } else {
    assert(from > to);
    n = (unsigned int) (from - to) + 1;
    by = -1;
  }
  
  if(n == 0) return NULL;
  
  s = (int*) malloc(sizeof(int) * n);
  assert(s);
  s[0] = from;
  for(i=1; i<n; i++) {
    s[i] = s[i-1] + by;
  }
  return s;
}

/*
 * allocates a new double array of size n1
 */
double* new_vector(unsigned int n)
{
  double *v;
  if(n == 0) return NULL;
  v = (double*) malloc(sizeof(double) * n);
  return v;
}

int* new_ivector(int n)
{
  int* v;
  v = (int*) malloc(sizeof(int)*n);
  return v;
}

int* new_zero_ivector(int n)
{
  int *v,i;
  v = new_ivector(n);
  for(i=0;i<n;i++)
    v[i] = 0;
  return v;
}

/*
 * allocates a new double array of size n1
 * and fills it with zeros
 */
double* new_zero_vector(unsigned int n)
{
  double *v;
  v = new_vector(n);
  zerov(v, n);
  return v;
}

/*
 * allocates a new double array of size n1
 * and fills it with the contents of vold
 */
double* new_dup_vector(double* vold, unsigned int n)
{
  double *v;
  v = new_vector(n);
  dupv(v, vold, n);
  return v;
}

void dupiv(int* v, int* vold, int n)
{
  int i;
  for(i=0;i<n;i++)
    v[i] = vold[i];
}

/*
 * copies vold to v 
 * (assumes v has already been allcocated)
 */
void dupv(double *v, double* vold, unsigned int n)
{
  unsigned int i;
  for(i=0; i<n; i++) v[i] = vold[i];
}

/*
 * swaps the pointer of v2 to v1, and vice-versa
 * (aviods copying via dupv)
 */
void swap_vector(double **v1, double **v2)
{
  double *temp;
  temp = *v1;
  *v1 = *v2;
  *v2 = temp;
}

/*
 * zeros out v
 * (assumes that it has already been allocated)
 */
void zerov(double*v, unsigned int n)
{
  unsigned int i;
  for(i=0; i<n; i++) v[i] = 0;
}

/*
 * printing a vector out to outfile
 */
void printVector(double *v, unsigned int n, FILE *outfile)
{
  unsigned int i;
  for(i=0; i<n-1; i++) fprintf(outfile, "%.6f ", v[i]);
  fprintf(outfile, "%.6f\n",v[i]);
}

void printiVector(int *v, int n, FILE *outfile)
{
  int i;
  for(i=0;i<n-1;i++) fprintf(outfile, "%d ", v[i]);
  fprintf(outfile, "%d\n", v[i]);
}

void printVectorR(double *v, int n, FILE *outfile)
{
  int i;
  fprintf(outfile,"c(");
  for(i=0; i<n-1; i++) 
    fprintf(outfile, "%.6f, ", v[i]);
  fprintf(outfile, "%.6f)\n",v[i]);
}

void printiVectorR(int *v, int n, FILE *outfile)
{
  int i;
  fprintf(outfile,"c(");
  for(i=0;i<n-1;i++) fprintf(outfile, "%d, ", v[i]);
  fprintf(outfile, "%d)\n", v[i]);
}

void printMatrixR(double** M, int d1, int d2, FILE *outfile)
{
  int i,j;
  fprintf(outfile,"matrix(c(");
  for(i=0;i<d1;i++)
    for(j=0;j<d2;j++)
      fprintf(outfile,"%.6f, ",M[i][j]);
  fprintf(outfile,"),nrow=%d,ncol=%d,byrow=TRUE)\n",d1,d2);
}

/*      
 * return the minimum element in the vector.
 * pass back the index of the minimum through 
 * the which pointer
 */
double min(double *v, unsigned int n, unsigned int *which)
{
  unsigned int i;
  double min;
  
  *which = 0;
  min = v[0];     
  
  for(i=1; i<n; i++) {    
    if(v[i] < min)  {
      min = v[i];
      *which = i;
    }
  }
  return min;
}

/*
 * return the maximum element in the vector.
 * pass back the index of the maximum through 
 * the which pointer
 */
double max(double *v, unsigned int n, unsigned int *which)
{
  unsigned int i;
  double max;
  
  *which = 0;
  max = v[0];
  
  for(i=1; i<n; i++) {
    if(v[i] > max)  {
      max = v[i];
      *which = i;
    }
  }
  
  return max;
}
