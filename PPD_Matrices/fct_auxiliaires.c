#include "header.h"

void remplis_colonne_matrice_avec_vecteur(gsl_matrix *matrice, gsl_vector *vecteur, size_t numero_colonne)
{
    if (matrice->size1 != vecteur->size)
    {
        printf("erreur dans remplis_1e_colonne_matrice_avec_vecteur. La matrice et le vecteur en entrée n'ont pas la même taille.\n");
    }
    else
    {
        for (size_t i = 0; i < vecteur->size; i++)
        {
            gsl_matrix_set(matrice, i, numero_colonne, vecteur->data[i]);
        }
    }
}

double produit_scalaire(gsl_vector *yk, gsl_vector *yk_suivant)
{
    double res = 0.0;

    for (size_t k = 0; k < yk->size; k++)
    {
        res += yk->data[k] * yk_suivant->data[k];
    }
    return res;
}

// effectue les produits scalaires et les stocke dans la matrice B1 et B0
// B0 correspond à Bm-1, B1 correspond à Bm dans l'énoncé
// écrase le yk par le yk+1 à chaque tour
void projection(gsl_spmatrix *A, gsl_matrix *B0, gsl_matrix *B1, gsl_matrix *Vm, gsl_vector *yk, size_t taille_sous_espace)
{
    printf("\ndimension matrice: %ld x %ld\n", A->size1, A->size2);

    double Ck = 0.0, temp = 0.0;
    int i, j, k, l, incremente_quand_k_pair = 0;
    gsl_vector *yk_suivant = gsl_vector_alloc(yk->size);
    printf("\ntaille sous espace: %ld\n\n", taille_sous_espace);

    Ck = produit_scalaire(yk, yk); // calcul de C0
    gsl_matrix_set(B0, 0, 0, Ck);  // stocke C0 dans B0
    remplis_colonne_matrice_avec_vecteur(Vm, yk, 0);

    // pour chaque index k de Bk
    for (k = 1; k <= 2 * taille_sous_espace - 1; k++)
    {
        if (k % 2 == 0) // si k est pair alors Ck=produit_scalaire(yk, yk)
        {
            Ck = produit_scalaire(yk, yk);
            incremente_quand_k_pair++;
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

                if (k < taille_sous_espace + incremente_quand_k_pair) // on stocke yk dans Vm
                {

                    remplis_colonne_matrice_avec_vecteur(Vm, yk_suivant, k - incremente_quand_k_pair);
                }
                temp = 0.0;
            }

            // calcul de Ck = produit scalaire(yk, yk_suivant)
            Ck = produit_scalaire(yk, yk_suivant);

            // copie les éléments de yk_suivant dans yk pour le prochain tour de boucle
            gsl_vector_memcpy(yk, yk_suivant);
        }

        i = k;

        if (i < taille_sous_espace)
        {
            gsl_matrix_set(B0, i, 0, Ck); // on remplit la première colonne de B0
        }

        if (k <= taille_sous_espace)
        {
            j = k;
            i = 0;
            while (j >= 1)
            {

                gsl_matrix_set(B1, i, j - 1, Ck); // on stocke Ck dans la matrice B1

                if (k != taille_sous_espace)
                {
                    gsl_matrix_set(B0, i, j, Ck); // on stocke Ck dans la matrice B0
                }
                else
                {
                    if (j != 1) // on ne remplit pas la première colonne de B0 car on la remplit déjà plus haut
                    {
                        gsl_matrix_set(B0, i + 1, j - 1, Ck); // on stocke Ck dans la matrice B0
                    }
                }

                i++;
                j--;
            }
        }
        else
        {
            j = taille_sous_espace;
            i = k - taille_sous_espace;
            while (j >= 1 + k - taille_sous_espace)
            { // on ne remplit pas la première colonne de B0 car on la remplit déjà plus haut
                if (i != taille_sous_espace - 1)
                { // on ne remplit pas la dernière colonne non plus d'où le j+1 et j-1
                    printf("fichier %s ligne %d\n", __FILE__, __LINE__);
                    printf("valeur de k: %d\n", k);
                    gsl_matrix_set(B0, i + 1, j - 1, Ck); // on stocke Ck dans la matrice B0
                    printf("fichier %s ligne %d\n", __FILE__, __LINE__);
                }
                gsl_matrix_set(B1, i, j - 1, Ck); // on stocke Ck dans la matrice B1
                j--;
                i++;
            }
        }
    }
    printf("B1:\n");
    affiche_matrice(B1);
    printf("B0:\n");
    affiche_matrice(B0);
    printf("Vm:\n");
    affiche_matrice(Vm);
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
    printf("\n\n");
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

void calcule_valeurs_propre(gsl_matrix *matrix, gsl_vector* valeurs_propres, gsl_matrix * vecteurs_propres)
{
    
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (matrix->size1);
    gsl_vector_complex *eval = gsl_vector_complex_alloc (matrix->size1);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc (matrix->size1, matrix->size2);
    gsl_eigen_nonsymmv (matrix, eval, evec, w);
    gsl_eigen_nonsymmv_free (w);


    for (int i = 0; i < matrix->size1; i++)
    {
        gsl_complex eval_i = gsl_vector_complex_get (eval, i);
        gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);

        printf ("Valeur propre = %g\n", GSL_REAL(eval_i));
        printf ("Vecteurs propres : \n");
        for (int j = 0; j < matrix->size1; ++j)
        {
            gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
            printf("%g\n", GSL_REAL(z));
        }
    }
    //Pour les valeurs propres 
    // traduit le gsl_vector_complex en gsl_vector
    for (size_t i = 0; i < eval->size; i++) {
        gsl_vector_set (valeurs_propres, i, GSL_REAL(gsl_vector_complex_get(eval, i)));
    }

    //Pour les vecteurs propres
    // traduit gsl_matrix_complex en gsl_matrix
    for (size_t i = 0; i < evec->size1; i++) {
        for (size_t j = 0; j < evec->size2; j++) {
            gsl_matrix_set(vecteurs_propres, i, j, GSL_REAL(gsl_matrix_complex_get(evec, i, j)));
        }
    }

    //affiche_matrice(vecteurs_propres);


  
    gsl_vector_complex_free (eval);
    gsl_matrix_complex_free (evec);
}


/*void calcule_vecteurs_propres(gsl_matrix *matrix, gsl_vector* valeurs_propres, gsl_matrix * vecteurs_propres)
{
    // Boucle pour calculer les vecteurs propres pour chaque valeur propre
    for (int i = 0; i < matrix->size1; i++) {
        double eigenvalue = gsl_vector_get(valeurs_propres, i);

        // Initialisation du vecteur propre
        gsl_vector* eigenvector = gsl_vector_alloc(matrix->size1);

        // Initialisation de la matrice A - λI
        gsl_matrix* A_lambdaI = gsl_matrix_alloc(matrix->size1, matrix->size2);
        gsl_matrix_memcpy(A_lambdaI, matrix);
        gsl_matrix_add_constant(A_lambdaI, -eigenvalue);
        gsl_permutation* p = gsl_permutation_alloc(matrix->size1);
        // Résolution du système linéaire (A - λI)x = 0
        int signum;
        
        // LU decomposition
        gsl_linalg_HH_decomp(A_lambdaI, p, &signum);
        int status = gsl_linalg_HH_solve(A_lambdaI, p, eigenvector, eigenvector);
        if (status) {
            printf("Erreur lors de la résolution du système linéaire\n");
        }

        // Stockage du vecteur propre dans la matrice des vecteurs propres
        gsl_matrix_set_col(vecteurs_propres, i, eigenvector);

        // Libération de la mémoire
        gsl_matrix_free(A_lambdaI);
        gsl_vector_free(eigenvector);
    }

    // Affichage des vecteurs propres
    for (int i = 0; i < matrix->size1; i++) {
        for (int j = 0; j < matrix->size2; j++) {
            printf("%g ", gsl_matrix_get(vecteurs_propres, i, j));
        }
        printf("\n");
    }

}
*/

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

gsl_matrix *multiplie_matrices(gsl_matrix *matrice1, gsl_matrix *matrice2)
{
    if (matrice1->size2 != matrice2->size1)
    {
        printf("le nombre de colonnes de matrice1 est différent du nombre de lignes de matrice2.La multiplication a échoué\n");
        return matrice1;
    }

    gsl_matrix *resultat = gsl_matrix_alloc(matrice1->size1, matrice2->size2);
    double temp = 0.0;
    for (int i = 0; i < matrice1->size1; i++)
    {
        for (int j = 0; j < matrice2->size2; j++)
        {
            for (int k = 0; k < matrice1->size2; k++)
            {
                temp += gsl_matrix_get(matrice1, i, k) * gsl_matrix_get(matrice2, k, j);
            }
            gsl_matrix_set(resultat, i, j, temp);
            temp = 0;
        }
    }
    return resultat;
}
