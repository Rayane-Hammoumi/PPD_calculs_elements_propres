#include "header.h"

// effectue les produits scalaires et les stocke dans la matrice B
// écrase le yk par le yk+1 à chaque tour
double produit_scalaire(gsl_vector *yk, gsl_vector *yk_suivant)
{
    double res = 0.0;

    for (size_t k = 0; k < yk->size; k++)
    {
        res += yk->data[k] * yk_suivant->data[k];
    }
    return res;
}

void projection(gsl_spmatrix *A, gsl_matrix *B, gsl_vector *yk, size_t taille_sous_espace)
{
    printf("\ndimension matrice: %ld x %ld\n", A->size1, A->size2);

    double Ck = 0.0, temp = 0.0;
    int i, j, k, l;
    gsl_vector *yk_suivant = gsl_vector_alloc(yk->size);
    printf("taille sous espace: %ld\n", taille_sous_espace);

    // pour chaque index k de Bk
    for (k = 1; k <= 2 * taille_sous_espace - 1; k++)
    {
        if (k % 2 == 0) // si k est pair alors Ck=produit_scalaire(yk, yk)
        {
            Ck = produit_scalaire(yk, yk);
            printf("Ck: %g\n", Ck);
        }
        else // sinon Ck=produit_scalaire(yk, yk_suivant)
        {
            // pour chaque élément l de yk_suivant
            for (l = 0; l < A->size1; l++)
            {
                // on calcule l'élément l
                for (i = 0; i < (A->size1); i++)
                {
                    temp += yk->data[i] * gsl_spmatrix_get(A, l, i);
                }
                yk_suivant->data[l] = temp;

                // printf("\nvaleur de temp: %g\n", temp);
                temp = 0.0;
            }

            // calcul de Ck = produit scalaire(yk, yk_suivant)
            Ck = produit_scalaire(yk, yk_suivant);
            printf("Ck: %g\n", Ck);

            // copie les éléments de yk_suivant dans yk pour le prochain tour de boucle
            gsl_vector_memcpy(yk, yk_suivant);
        }

        if (k <= taille_sous_espace)
        {
            j = k;
            i = 0;
            while (j >= 1)
            {
                printf("%d\n", __LINE__);
                gsl_matrix_set(B, i, j - 1, Ck); // on stocke Ck dans la matrice B
                i++;
                j--;
            }
        }
        else
        {
            j = taille_sous_espace;
            i = k - taille_sous_espace;
            while (j >= 1 + k - taille_sous_espace)
            {
                printf("%d\n", __LINE__);
                gsl_matrix_set(B, i, j - 1, Ck); // on stocke Ck dans la matrice B
                j--;
                i++;
            }
        }
    }
    affiche_matrice(B);
    gsl_vector_free(yk_suivant);
}

void affiche_matrice(gsl_matrix *A)
{
    for (int i = 0; i < A->size1; i++)
    {
        for (int j = 0; j < A->size2; j++)
        {
            printf(" %g", gsl_matrix_get(A, i, j));
        }
        printf("\n");
    }
}

gsl_spmatrix *lit_fichier_mat(char nomFichier[])
{
    FILE *fp;
    fp = fopen(nomFichier, "r");

    gsl_spmatrix *A = gsl_spmatrix_fscanf(fp);

    fclose(fp);

    // affiche la matrice importée
    // gsl_spmatrix_fprintf(stdout, A, "%.g");

    return A;
}

void calcule_valeurs_propre(gsl_matrix *matrix)
{
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(matrix->size1);
    gsl_vector *eval = gsl_vector_alloc(matrix->size1);
    gsl_matrix *evec = gsl_matrix_alloc(matrix->size1, matrix->size1);
    gsl_eigen_symmv(matrix, eval, evec, w);
    gsl_eigen_symmv_free(w);

    // afficher les valeurs propres
    for (int i = 0; i < matrix->size1; i++)
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
