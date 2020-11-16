#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/* Prototypes for the Fortran functions */

/* Dot product */
double dot_blas_(const int *n, const double *x1, const double *x2);
double dot_for_(const int *n, const double *x1, const double *x2);
double dot_intrinsic_(const int *n, const double *x1, const double *x2);

/* Matrix - vector multiply */
void matvec_blas_(const int *m, const int *n, const double *mat, const double *vec, double *res);
void matvec_for_(const int *m, const int *n, const double *mat, const double *vec, double *res);
void matvec_intrinsic_(const int *m, const int *n, const double *mat, const double *vec, double *res);

/* Matrix - matrix multiply */
void matmat_blas_(const int *m, const int *n, const int *k, const double *a, const double *b, double *c);
void matmat_for_(const int *m, const int *n, const int *k, const double *a, const double *b, double *c);
void matmat_intrinsic_(const int *m, const int *n, const int *k, const double *a, const double *b, double *c);

double time_elapsed(const struct timespec start, const struct timespec end) {
  return (end.tv_sec-start.tv_sec)+1.0e-9*(end.tv_nsec-start.tv_nsec);
}

int main(void) {
  /* Size of array and vector index variable */
  int N, i, j;
  /* Timing stuff */
  struct timespec start, stop;
  /* Vector, and dot products */
  double *x=NULL, dblas, dfor, dint;

  /* Matrix, vector and result */
  double *mat=NULL, *vec=NULL, *res=NULL;

  /* Matrix, matrix, matrix */
  double *a=NULL, *b=NULL, *c=NULL;

  /* Times of computing dot products */
  double tblas, tfor, tint;

  printf("Dot products:\n");
  //  for(N=10;N<=1e9;N*=10) {
  for(N=10;N<=1e3;N*=10) {
    /* Increase size of x */
    x=realloc(x,N*sizeof(double));
    /* Initialize x */
    for(i=0;i<N;i++)
      x[i]=i*0.1;

    /* Compute dot products */

    /* BLAS */
    clock_gettime(CLOCK_REALTIME,&start);
    dblas=dot_blas_(&N,x,x);
    clock_gettime(CLOCK_REALTIME,&stop);
    tblas=time_elapsed(start,stop);

    /* For loop */
    clock_gettime(CLOCK_REALTIME,&start);
    dfor=dot_for_(&N,x,x);
    clock_gettime(CLOCK_REALTIME,&stop);
    tfor=time_elapsed(start,stop);

    /* Intrinsic Fortran */
    clock_gettime(CLOCK_REALTIME,&start);
    dint=dot_intrinsic_(&N,x,x);
    clock_gettime(CLOCK_REALTIME,&stop);
    tint=time_elapsed(start,stop);

    /* Print out results */
    printf("N = %8i\n",N);
    printf("\tblas\t%e\t%e\n",dblas,tblas);
    printf("\tfor\t%e\t%e\n",dfor,tfor);
    printf("\tintr\t%e\t%e\n",dint,tint);
  }
  free(x);

  /* Same thing for matrix-vector product */
  printf("\n\nMatrix-vector:\n");
  for(N=10;N<=1e4;N*=10) {
    //for(N=10;N<=1e2;N*=10) {
    /* Increase size of mat and vec */
    mat=realloc(mat,N*N*sizeof(double));
    vec=realloc(vec,N*sizeof(double));
    res=realloc(res,N*sizeof(double));

    /* Initialize mat and vec */
    for(j=0;j<N;j++)
      for(i=0;i<N;i++)
	mat[j*N+i]=1.0/(i+j+2);
    
    for(i=0;i<N;i++)
      vec[i]=i+1.0;

    /* Compute products */

    /* BLAS */
    clock_gettime(CLOCK_REALTIME,&start);
    matvec_blas_(&N,&N,mat,vec,res);
    clock_gettime(CLOCK_REALTIME,&stop);
    tblas=time_elapsed(start,stop);

    /* For loop */
    clock_gettime(CLOCK_REALTIME,&start);
    matvec_for_(&N,&N,mat,vec,res);
    clock_gettime(CLOCK_REALTIME,&stop);
    tfor=time_elapsed(start,stop);

    /* Intrinsic Fortran */
    clock_gettime(CLOCK_REALTIME,&start);
    matvec_intrinsic_(&N,&N,mat,vec,res);
    clock_gettime(CLOCK_REALTIME,&stop);
    tint=time_elapsed(start,stop);

    /* Print out results */
    printf("N = %8i\n",N);
    printf("\tblas\t%e\n",tblas);
    printf("\tfor\t%e\n",tfor);
    printf("\tintr\t%e\n",tint);
  }
  free(mat);
  free(vec);
  free(res);


  /* Same thing for matrix-matrix product */
  printf("\n\nMatrix-matrix:\n");
  for(N=10;N<=1e2;N*=10) {
    /* Increase size of variables */
    a=realloc(a,N*N*sizeof(double));
    b=realloc(b,N*N*sizeof(double));
    c=realloc(c,N*N*sizeof(double));

    /* Initialize */
    for(j=0;j<N;j++)
      for(i=0;i<N;i++) {
	a[j*N+i]=1.0/(i+j+2);
	b[j*N+i]=(i+j+2.0);
      }
    
    /* Compute products */

    /* BLAS */
    clock_gettime(CLOCK_REALTIME,&start);
    matmat_blas_(&N,&N,&N,a,b,c);
    clock_gettime(CLOCK_REALTIME,&stop);
    tblas=time_elapsed(start,stop);

    /* For loop */
    clock_gettime(CLOCK_REALTIME,&start);
    matmat_for_(&N,&N,&N,a,b,c);
    clock_gettime(CLOCK_REALTIME,&stop);
    tfor=time_elapsed(start,stop);

    /* Intrinsic Fortran */
    clock_gettime(CLOCK_REALTIME,&start);
    matmat_intrinsic_(&N,&N,&N,a,b,c);
    clock_gettime(CLOCK_REALTIME,&stop);
    tint=time_elapsed(start,stop);

    /* Print out results */
    printf("N = %8i\n",N);
    printf("\tblas\t%e\n",tblas);
    printf("\tfor\t%e\n",tfor);
    printf("\tintr\t%e\n",tint);
  }
  free(a);
  free(b);
  free(c);

  

  return 0;
}
