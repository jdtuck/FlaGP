#ifndef __COVAR_H__
#define __COVAR_H__

// changed by grant
void covar_symm(const int col, double **X, const int n,
                double d, double g, double **K);
void covar_symm_R(int *col_in, double *X_in, int *n_in,
                double *d_in, double *g_in, double *K_in);
void covar(const int col, double **X1, const int n1, double **X2,
	   const int n2, double d, double **K);
void covar_R(int *col_in, double *X1_in, int *n1_in, double *X2_in,
             int *n2_in, double *d_in, double *K_in);
//
#endif

