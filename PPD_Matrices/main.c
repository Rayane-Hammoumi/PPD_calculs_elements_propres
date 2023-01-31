#include "header.h"

int main(int argc, char *argv[])
{
  if (argc < 4) // si il n'y a pas assez d'arguments dans la commande
  {
    printf("Veuillez spécifier tous les arguments.\nUtilisation: ./calcul_elements_propres chemin_fichier taille_sous_espace précision\n");
    return 1;
  }

  else // s'il y a assez d'arguments dans la commande
  {
    char *verifie_strtol; // pointeur utilisé par strtol pour vérifier que la taille spécifiée est bien un nombre
    size_t taille_sous_espace = strtol(argv[2], &verifie_strtol, 10);
    if (*verifie_strtol != '\0')
    {
      printf("La taille spécifiée pour le sous-espace n'est pas un nombre\nUtilisation: ./calcul_elements_propres chemin_fichier taille_sous_espace précision\n");
      return 1;
    }

    double epsilon;
    sscanf(argv[3], "%lf", &epsilon);

    // déclarations et allocations de mémoire
    gsl_spmatrix *A = lit_fichier_mat(argv[1]); // on importe la matrice du fichier
    gsl_matrix *B0 = gsl_matrix_calloc(taille_sous_espace, taille_sous_espace);
    gsl_matrix *B1 = gsl_matrix_calloc(taille_sous_espace, taille_sous_espace);
    gsl_matrix *Vm = gsl_matrix_calloc(A->size1, taille_sous_espace);
    gsl_vector *yk = gsl_vector_calloc(A->size1);
    gsl_vector *valeurs_propres = gsl_vector_calloc(taille_sous_espace);                     // on stocke les valeurs propres dans un gsl_vector
    gsl_matrix *vecteurs_propres = gsl_matrix_alloc(taille_sous_espace, taille_sous_espace); // on stocke les vecteurs propres dans une gsl_matrix
    gsl_matrix *qi = gsl_matrix_alloc(A->size1, taille_sous_espace);                         // matrice qui contient les vecteurs qi
    gsl_vector *result = gsl_vector_alloc(A->size1);

    int precision_atteinte = 0, compteur_iterations = 0;

    // pour stocker le temps d'exécution du code
    double start_time, end_time, time;

    gsl_vector_set_zero(yk);  // le y0 de la première itération est égal à la base (1, 0, 0...)
    gsl_vector_set(yk, 0, 1); // en effet on choisit de prendre x = (1, 0, 0...) donc ||x||=1

    start_time = omp_get_wtime();

    while (!precision_atteinte)
    {
      double start_timepr, end_timepr, timepr;

      start_timepr = omp_get_wtime();

      printf("\nProjection :\n");
      // calcul de Bm (<=>B1), de Bm-1 (<=>B0) et de Vm
      projection(A, B0, B1, Vm, yk, taille_sous_espace); // remplit B0, B1 et Vm et les affiche

      end_timepr = omp_get_wtime();

      timepr = end_timepr - start_timepr;

      // calcul et affichage de Em
      gsl_matrix *Em = inverse_matrix(B0);
      printf("EM:\n");
      affiche_matrice(Em);

      // calcul et affichage de Fm
      gsl_matrix *Fm = multiplie_matrices(Em, B1);
      printf("Fm:\n");
      affiche_matrice(Fm);

      // calcul des valeurs/vecteurs propres de Fm
      printf("matrice vecteurs_propres:\n");
      calcule_valeurs_et_vecteurs_propre(Fm, valeurs_propres, vecteurs_propres);

      printf("vecteurs_propres size1 = %ld et Vm size1 = %ld\n", vecteurs_propres->size1, Vm->size1);
      printf("vecteurs_propres size2 = %ld et Vm size2 = %ld\n", vecteurs_propres->size2, Vm->size2);

      // Partie 3
      printf("\nProduit matrice-vecteur :\n");
      // pour stocker le temps d'exécution du code
      double start_time, end_time, time;

      start_time = omp_get_wtime();
      calcule_qi(A, qi, vecteurs_propres, Vm, taille_sous_espace);
      end_time = omp_get_wtime();

      time = end_time - start_time;

      affiche_matrice(qi);

      if (verifie_si_precision_atteinte(A, valeurs_propres, qi, epsilon) == 1)
      {
        precision_atteinte = 1;
      }
      else
      {
        gsl_matrix_get_col(yk, qi, qi->size2 - 1);                // on prend x = le dernier vecteur propre dans qi (sa dernière colonne)
        produit_constante_vecteur(1 / calcule_norme(yk), yk, yk); // yk devient x/||x||
        gsl_matrix_free(Em);
        gsl_matrix_free(Fm);
        compteur_iterations++;
        system("clear");
      }

      printf("----[Derniere iteration]----\n");
      if (timepr < 1)
        printf("[Projection] Temps d'execution : %f ms\n", (timepr * 1000.0));
      else
        printf("[Projection] Temps d'execution : %f s\n", timepr);

      if (time < 1)
        printf("[Produit matrice-vecteur] Temps d'execution : %f ms\n", (time * 1000.0));
      else
        printf("[Produit matrice-vecteur] Temps d'execution : %f s\n", time);
    }
    printf("\n");
    //  libération de la mémoire allouée
    gsl_spmatrix_free(A);
    gsl_matrix_free(B0);
    gsl_matrix_free(B1);
    gsl_matrix_free(Vm);
    gsl_vector_free(yk);
    gsl_vector_free(valeurs_propres);
    gsl_matrix_free(vecteurs_propres);
    gsl_matrix_free(qi);
    gsl_vector_free(result);

    end_time = omp_get_wtime();
    time = end_time - start_time;
    if (time < 1)
    {
      printf("[Temps d'execution TOTAL] = %f ms\n", (time * 1000.0));
    }
    else
    {
      printf("[Temps d'execution TOTAL] = %f s\n", time);
      printf("Nombre d'itérations: %d\n", compteur_iterations);
      printf("Precision: %g\n", epsilon);
      printf("Nombre de Threads: %d\n", omp_get_max_threads());
    }
  }
  /*
    double data[] = { 0, 3, 5,
                      5, 5, 2,};

    gsl_matrix_view m = gsl_matrix_view_array (data, 2, 3);
    affiche_matrice(&m.matrix);
    gsl_vector* v = gsl_vector_calloc(m.matrix.size2);
    gsl_vector_set(v, 0, 3);
    gsl_vector_set(v, 1, 4);
    gsl_vector_set(v, 2, 3);

    produit_matrice_vecteur(&m.matrix, v);*/

  return 0;
}