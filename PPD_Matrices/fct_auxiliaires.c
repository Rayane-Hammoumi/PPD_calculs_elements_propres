#include "header.h"

// renvoie 1 si la précision est atteinte, sinon 0
int verifie_si_precision_atteinte(gsl_spmatrix *A, gsl_vector *valeurs_propres, gsl_matrix *qi, double epsilon)
{
    int max = 0;
    double temp = 0.0;
    gsl_vector *A_multiplie_par_qi = gsl_vector_alloc(A->size1);
    gsl_vector *lambda_multiplie_par_qi = gsl_vector_alloc(qi->size1);
    gsl_vector *res_soustraction = gsl_vector_alloc(qi->size1);

    //parallélisation de la boucle for
    #pragma omp parallel for
    for (int i = 0; i < qi->size2; i++)
    {
        gsl_vector_view tempui = gsl_matrix_column(qi, i); // on récupère une colonne dans la matrice contenant les vecteurs qi
        gsl_vector *vecteur_qi = &tempui.vector;

        produit_spmatrice_vecteur(A, vecteur_qi, A_multiplie_par_qi);
        produit_constante_vecteur(gsl_vector_get(valeurs_propres, i), vecteur_qi, lambda_multiplie_par_qi);
        soustrait_vecteur2_au_vecteur1(A_multiplie_par_qi, lambda_multiplie_par_qi, res_soustraction);
        temp = calcule_norme(res_soustraction);

        if (temp > max)
        {
            max = temp;
        }
    }

    if (max > epsilon)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}


void soustrait_vecteur2_au_vecteur1(gsl_vector *vecteur1, gsl_vector *vecteur2, gsl_vector *resultat)
{
    if (vecteur1->size != vecteur2->size)
    {
        printf("Échec de la soustraction des vecteurs. Ils n'ont pas la même taille.\n");
        exit(1);
    }

    //parallélisation de la boucle for
    #pragma omp parallel for
    for (int i = 0; i < vecteur1->size; i++)
    {
        gsl_vector_set(resultat, i, gsl_vector_get(vecteur1, i) - gsl_vector_get(vecteur2, i));
    }
}


void produit_constante_vecteur(double constante, gsl_vector *vecteur, gsl_vector *resultat)
{
    //parallélisation de la boucle for
    #pragma omp parallel for
    for (int i = 0; i < vecteur->size; i++)
    {
        gsl_vector_set(resultat, i, gsl_vector_get(vecteur, i) * constante);
    }
}


double calcule_norme(gsl_vector *vecteur)
{
    double resultat = 0;

    //parallélisation de la boucle for
    #pragma omp parallel for reduction(+:resultat)
    for (int i = 0; i < vecteur->size; i++)
    {
        resultat += vecteur->data[i] * vecteur->data[i];
    }

    resultat = sqrt(resultat);

    return resultat;
}


double produit_scalaire(gsl_vector *yk, gsl_vector *yk_suivant)
{
    double res = 0.0;

    double start_time, end_time, time;
 
    start_time = omp_get_wtime();


    // pour chaque élément de result
    #pragma omp parallel for reduction(+:res)
    for (size_t k = 0; k < yk->size; k++)
    {
        res += yk->data[k] * yk_suivant->data[k];
    }

    end_time = omp_get_wtime();
    time = end_time - start_time;
    if(time < 1)
      printf("[produit_scalaire] Temps d'execution : %f ms\n", (time*1000.0));
    else
      printf("[produit_scalaire] Temps d'execution : %f s\n", time);


    return res;
}
// effectue les produits scalaires et les stocke dans la matrice B1 et B0
// B0 correspond à Bm-1, B1 correspond à Bm dans l'énoncé
// écrase le yk par le yk+1 à chaque tour
// remplit la matrice Vm
void projection(gsl_spmatrix *A, gsl_matrix *B0, gsl_matrix *B1, gsl_matrix *Vm, gsl_vector *yk, size_t taille_sous_espace)
{
    printf("\ndimension matrice: %ld x %ld\n", A->size1, A->size2);

    double Ck = 0.0;
    int i, j, k, incremente_quand_k_pair = 0;
    gsl_vector *yk_suivant = gsl_vector_alloc(yk->size);
    printf("\ntaille sous espace: %ld\n\n", taille_sous_espace);

    Ck = produit_scalaire(yk, yk); // calcul de C0
    gsl_matrix_set(B0, 0, 0, Ck);  // stocke C0 dans B0
    gsl_matrix_set_col(Vm, 0, yk);

    // pour chaque index k de Bk
    #pragma omp for schedule(static)
    for (k = 1; k <= 2 * taille_sous_espace - 1; k++)
    {
        if (k % 2 == 0) // si k est pair alors Ck=produit_scalaire(yk, yk)
        {
            Ck = produit_scalaire(yk, yk);
            incremente_quand_k_pair++;
        }
        else // sinon Ck=produit_scalaire(yk, yk_suivant)
        {
            produit_spmatrice_vecteur(A, yk, yk_suivant);

            if (k < taille_sous_espace + incremente_quand_k_pair) // on stocke yk dans Vm
            {
                gsl_matrix_set_col(Vm, k - incremente_quand_k_pair, yk_suivant);
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
void affiche_vecteur(gsl_vector *v, int taille)
{
    // printf("Vecteur :\n");
    for (int i = 0; i < taille; i++)
    {
        printf("%g ", gsl_vector_get(v, i));
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

void calcule_valeurs_et_vecteurs_propre(gsl_matrix *matrix, gsl_vector *valeurs_propres, gsl_matrix *vecteurs_propres)
{

    gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(matrix->size1);
    gsl_vector_complex *eval = gsl_vector_complex_alloc(matrix->size1);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(matrix->size1, matrix->size2);
    gsl_eigen_nonsymmv(matrix, eval, evec, w);
    gsl_eigen_nonsymmv_free(w);

    for (int i = 0; i < matrix->size1; i++)
    {
        gsl_complex eval_i = gsl_vector_complex_get(eval, i);
        gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec, i);

        printf("Valeur propre = %g\n", GSL_REAL(eval_i));
        printf("Vecteurs propres : \n");
        for (int j = 0; j < matrix->size1; ++j)
        {
            gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
            printf("%g\n", GSL_REAL(z));
        }
    }
    // Pour les valeurs propres
    //  traduit le gsl_vector_complex en gsl_vector
    for (size_t i = 0; i < eval->size; i++)
    {
        gsl_vector_set(valeurs_propres, i, GSL_REAL(gsl_vector_complex_get(eval, i)));
    }

    // Pour les vecteurs propres
    //  traduit gsl_matrix_complex en gsl_matrix
    //  on met les vecteurs propres associés a une valeur propre ligne par ligne
    for (size_t i = 0; i < evec->size1; i++)
    {
        for (size_t j = 0; j < evec->size2; j++)
        {
            gsl_matrix_set(vecteurs_propres, i, j, GSL_REAL(gsl_matrix_complex_get(evec, i, j)));
        }
    }

    affiche_matrice(vecteurs_propres);

    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
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

void produit_matrice_vecteur(gsl_matrix *m, gsl_vector *v, gsl_vector *resultat)
{
    // printf("Produit matrice-vecteur :\n");

    double temp = 0.0;

    // pour chaque élément de resultat
    #pragma omp for schedule(static)
    for (int i = 0; i < m->size1; i++)
    {
        // on le calcule (somme des produits)
        for (int j = 0; j < m->size2; j++)
        {
            temp += gsl_matrix_get(m, i, j) * gsl_vector_get(v, j);
        }
        gsl_vector_set(resultat, i, temp);
        temp = 0.0;
    }

    // affiche_vecteur(resultat, m->size1);
}

void produit_spmatrice_vecteur(gsl_spmatrix *m, gsl_vector *v, gsl_vector *resultat)
{
    // printf("Produit matrice-vecteur :\n");

    double temp = 0.0;

    // pour chaque élément de resultat
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < m->size1; i++)
    {
        // on le calcule (somme des produits)
        for (int j = 0; j < m->size2; j++)
        {
            temp += gsl_spmatrix_get(m, i, j) * gsl_vector_get(v, j);
        }
        gsl_vector_set(resultat, i, temp);
        temp = 0.0;
    }

    // affiche_vecteur(result, m->size1);
}

void calcule_qi(gsl_spmatrix *A, gsl_matrix *qi, gsl_matrix *vecteurs_propres, gsl_matrix *Vm, size_t taille_sous_espace)
{

    gsl_vector *result = gsl_vector_alloc(A->size1);

    // calculs des vecteurs qi qu'on stocke dans la gsl_matrix qi
    #pragma omp for schedule(static)
    for (size_t i = 0; i < taille_sous_espace; i++)
    {
      gsl_vector_view tempui = gsl_matrix_column(vecteurs_propres, i); // vecteurs propres d'une valeur propre (tt la colonne)
      gsl_vector *ui = &tempui.vector;
      // affiche_vecteur(ui, taille_sous_espace);
      produit_matrice_vecteur(Vm, ui, result);
      gsl_matrix_set_col(qi, i, result);

    }
    
}