/****************************************************************************
 *
 * Local Approximate Gaussian Process Regression
 * Copyright (C) 2013, The University of Chicago
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 *
 * Questions? Contact Robert B. Gramacy (rbg@vt.edu)
 *
 ****************************************************************************/


#include "matrix.h"
#include <stdlib.h>
#ifdef RPRINT
#include <R.h>
#else
#include <math.h>
#endif


/*
 * covar_symm_R:
 *
 * calculate the correlation (K) between X1 and X2 with
 * an isotropic power exponential correlation function
 * with range d and nugget g -- assumes symmetric matrix
 */

void covar_symm(const int col, double **X, const int n,
                double d, double g, double **K)
{
  int i, j, k;

  /* calculate the covariance */
  for(i=0; i<n; i++) {
    for(j=i+1; j<n; j++) {
      K[i][j] = 0.0;
      for(k=0; k<col; k++) K[i][j] += sq(X[i][k] - X[j][k]);
      K[j][i] = K[i][j] = exp(0.0 - K[i][j]/d);
    }
    K[i][i] = 1.0 + g;
  }
}

void covar_symm_R(int *col_in, double *X_in, int *n_in,
                double *d_in, double *g_in, double *K_in)
{
  double **X;
  double **K;

  X = new_matrix_bones(X_in, *n_in, *col_in);
  K = new_matrix_bones(K_in, *n_in, *n_in);
  covar_symm(*col_in, X, *n_in, *d_in, *g_in, K);
  //double **Kchol = new_dup_matrix(K, *n_in, *n_in);
  //int info = linalg_dposv(*n_in, Kchol, K);
  //if(info) {
//#ifdef UNDEBUG
  //  printMatrix(K, *n_in, *n_in, stdout);
//#endif
  //  error("bad Cholesky decomp (info=%d)",
  //        info);
  //}
  free(X);
  free(K);
  //free(Kchol);
}


/*
 * covar:
 *
 * calculate the correlation (K) between X1 and X2 with
 * an isotropic power exponential correlation function
 * with range d and nugget g
 */
void covar(const int col, double **X1, const int n1, double **X2,
	   const int n2, double d, double **K)
{
  int i, j, k;

  /* calculate the covariance */
  for(i=0; i<n1; i++)
    for(j=0; j<n2; j++) {
      K[i][j] = 0.0;
      for(k=0; k<col; k++) K[i][j] += sq(X1[i][k] - X2[j][k]);
      K[i][j] = exp(0.0 - K[i][j]/d);
    }
}
void covar_R(int *col_in, double *X1_in, int *n1_in, double *X2_in,
             int *n2_in, double *d_in, double *K_in)
{
  double **X1;
  double **X2;
  double **K;

  X1 = new_matrix_bones(X1_in, *n1_in, *col_in);
  X2 = new_matrix_bones(X2_in, *n2_in, *col_in);
  K = new_matrix_bones(K_in, *n1_in, *n2_in);
  covar(*col_in, X1, *n1_in, X2, *n2_in, *d_in, K);
  free(X1);
  free(X2);
  free(K);
}
