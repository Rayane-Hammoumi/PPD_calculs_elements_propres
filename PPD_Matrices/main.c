#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

int main()
{
  int n = 3;
  gsl_matrix *A = gsl_matrix_alloc(n, n);

  /* on remplit la matrice A avec des valeurs */
  gsl_matrix_set(A, 0, 0, 2.0);
  gsl_matrix_set(A, 0, 1, 1.0);
  gsl_matrix_set(A, 0, 2, 0.0);
  gsl_matrix_set(A, 1, 0, 1.0);
  gsl_matrix_set(A, 1, 1, 2.0);
  gsl_matrix_set(A, 1, 2, 1.0);
  gsl_matrix_set(A, 2, 0, 0.0);
  gsl_matrix_set(A, 2, 1, 1.0);
  gsl_matrix_set(A, 2, 2, 2.0);

  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
  gsl_vector *eval = gsl_vector_alloc(n);
  gsl_matrix *evec = gsl_matrix_alloc(n, n);
  gsl_eigen_symmv(A, eval, evec, w);
  gsl_eigen_symmv_free(w);

  /* afficher les valeurs propres */
  for (int i = 0; i < n; i++)
  {
    printf("valeur propre = %g\n", gsl_vector_get(eval, i));
  }

  gsl_vector_free(eval);
  gsl_matrix_free(evec);
  gsl_matrix_free(A);
  return 0;
}
