#include "header.h"

int main(int argc, char *argv[])
{
  /*int n = 3;
  gsl_matrix *A = gsl_matrix_alloc(n, n);
  // on remplit la matrice A avec des valeurs

  gsl_matrix_set(A, 0, 0, 2.0);
  gsl_matrix_set(A, 0, 1, 1.0);
  gsl_matrix_set(A, 0, 2, 0.0);
  gsl_matrix_set(A, 1, 0, 1.0);
  gsl_matrix_set(A, 1, 1, 2.0);
  gsl_matrix_set(A, 1, 2, 1.0);
  gsl_matrix_set(A, 2, 0, 0.0);
  gsl_matrix_set(A, 2, 1, 1.0);
  gsl_matrix_set(A, 2, 2, 2.0);*/

  if (argc > 1)
  {

    gsl_spmatrix *sparseMatrix = lit_fichier_mat(argv[1]);

    // initialiser un objet gsl_matrix
    gsl_matrix *matrix = gsl_matrix_calloc(sparseMatrix->size1, sparseMatrix->size2);

    // parcourir tous les éléments de la matrice creuse
    for (int i = 0; i < sparseMatrix->size1; i++)
    {
      for (int j = sparseMatrix->i[i]; j < sparseMatrix->i[i + 1]; j++)
      {
        gsl_matrix_set(matrix, i, sparseMatrix->p[j], sparseMatrix->data[j]);
      }
    }

    int n = matrix->size1;

    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
    gsl_vector *eval = gsl_vector_alloc(n);
    gsl_matrix *evec = gsl_matrix_alloc(n, n);
    // gsl_eigen_symmv(A, eval, evec, w);
    gsl_eigen_symmv(matrix, eval, evec, w);
    gsl_eigen_symmv_free(w);

    // afficher les valeurs propres
    for (int i = 0; i < n; i++)
    {
      printf("valeur propre = %g\n", gsl_vector_get(eval, i));
    }

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    // gsl_matrix_free(A);
    gsl_matrix_free(matrix);
    gsl_spmatrix_free(sparseMatrix);
  }

  return 0;
}