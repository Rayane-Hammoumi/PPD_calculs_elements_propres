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
  gsl_matrix_set(A, 2, 2, 2.0);
  gsl_matrix_free(A);*/

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
    calcule_valeurs_propre(matrix);
    inverse_matrix(matrix);

    // Libérer la mémoire
    gsl_matrix_free(matrix);
    gsl_spmatrix_free(sparseMatrix);
  }

  return 0;
}