#include "header.h"

void affiche_matrice(gsl_matrix *A)
{
    for (int i = 0; i < A->size1; i++)
    {
        for (int j = 0; j < A->size2; j++)
        {
            printf("%g ", gsl_matrix_get(A, i, j));
        }
        printf("\n");
    }
}

gsl_spmatrix *lit_fichier_mat(char nomFichier[])
{
    FILE *fp;
    gsl_spmatrix *A;

    fp = fopen(nomFichier, "r");
    A = gsl_spmatrix_fscanf(fp);
    fclose(fp);

    // Do something with the matrix, e.g. print it
    gsl_spmatrix_fprintf(stdout, A, "%.g");

    return A;
}

void calcule_valeurs_propre(gsl_matrix *matrix)
{
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
}

gsl_matrix *inverse_matrix(gsl_matrix *A)
{
    int n = A->size1;
    // Allouer de la mémoire pour la matrice inverse
    gsl_matrix *inverse = gsl_matrix_alloc(n, n);

    // Allouer de la mémoire pour les vecteurs de permutation
    gsl_permutation *p = gsl_permutation_alloc(n);

    // Calculer l'inverse de la matrice
    int s;
    gsl_linalg_LU_decomp(A, p, &s);
    gsl_linalg_LU_invert(A, p, inverse);

    // Afficher la matrice inverse
    // affiche_matrice(inverse);

    // liberer la mémoire
    gsl_permutation_free(p);

    return inverse;
}
