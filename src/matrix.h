#ifdef __cplusplus
extern "C" {
#endif

#ifndef __MATRIX_H__
#define __MATRIX_H__ 

#include <stdio.h>

  typedef enum FIND_OP {LT=101, LEQ=102, EQ=103, GEQ=104, GT=105, NE=106} FIND_OP;
  
  /* my stuff */
  int ivequal(int* v1, int* v2, int n);
  int dvequal(double* v1, double* v2, int n);

  double** diag(double d, int n);
  
  void extractSubMatrix(double** A,int dim,int* ind,double** Asub);
  void insertSubMatrix(double** Asub, int dim, int* ind, double** A);

  int** new_2diarray(int d1, int d2);
  void rm_2diarray(int** x);
  double** new_2ddarray(int d1, int d2);
  void rm_2ddarray(double** x);
  float** new_2dfptrs(int d1, int d2, float* x);
  void rm_2dfptrs(float** x);
  int** new_2diptrs(int d1, int d2, int* x);
  void rm_2diptrs(int** x);
  double** new_2ddptrs(int d1, int d2, double* x);
  void rm_2ddptrs(double** x);
  float*** new_3dfarray(int  d1, int d2, int d3);
  void rm_3dfarray(float*** myarray, int d1);
  int*** new_3diarray(int  d1, int d2, int d3);
  void rm_3diarray(int*** myarray, int d1);
  double*** new_3ddarray(int  d1, int d2, int d3);
  void rm_3ddarray(double*** myarray, int d1);
  float*** new_3dfptrs(int  d1, int d2, int d3, float* x);
  void rm_3dfptrs(float*** myarray, int d1);
  int*** new_3diptrs(int  d1, int d2, int d3, int* x);
  void rm_3diptrs(int*** myarray, int d1);
  double*** new_3ddptrs(int  d1, int d2, int d3, double* x);
  void rm_3ddptrs(double*** myarray, int d1);
  double**** new_4ddptrs(int d1, int d2, int d3, int d4, double* x);
  void rm_4ddptrs(double**** x, int d1, int d2);
  double**** new_4ddarray(int d1, int d2, int d3, int d4);
  void rm_4ddarray(double**** x, int d1, int d2);
  void swap_3ddarray(double*** a, double*** b, int d1, int d2, int d3);

  
  /* bobby's stuff */
  void zero(double **M, unsigned int n1, unsigned int n2);
  void id(double **M, unsigned int n);
  double ** new_id_matrix(unsigned int n);
  double ** new_zero_matrix(unsigned int n1, unsigned int n2);
  double ** new_matrix(unsigned int m, unsigned int n);
  double ** new_dup_matrix(double** M, unsigned int n1, unsigned int n2);
  void dup_matrix(double** M1, double **M2, unsigned int n1, unsigned int n2);
  void swap_matrix(double **M1, double **M2, unsigned int n1, unsigned int n2);

  void delete_matrix(double** m);
  void printMatrix(double **M, unsigned int n, unsigned int col, FILE *outfile);
  double* ones(unsigned int n, double scale);
  double* dseq(double from, double to, double by);
  int* iseq(double from, double to);

  double* new_dup_vector(double* vold, unsigned int n);
  double* new_zero_vector(unsigned int n);
  double* new_vector(unsigned int n);

  int* new_ivector(int n);
  int* new_zero_ivector(int n);

  void printVector(double *v, unsigned int n, FILE *outfile);
  void printiVector(int* v, int n, FILE *outfile);
  void printVectorR(double* v, int n, FILE *outfile);
  void printiVectorR(int* v, int n, FILE *outfile);

  void dupv(double*v, double* vold, unsigned int n);
  void dupiv(int* v, int* vold, int n);

  void swap_vector(double **v1, double **v2);
  void zerov(double*v, unsigned int n);

  void printMatrixR(double** M, int d1, int d2, FILE *outfile);

  double max(double *v, unsigned int n, unsigned int *which);
  double min(double *v, unsigned int n, unsigned int *which);
  
#endif

#ifdef __cplusplus
}
#endif
