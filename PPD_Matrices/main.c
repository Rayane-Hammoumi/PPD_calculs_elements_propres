#include "header.h"

int main(int argc, char *argv[])
{
  if (argc < 3) // si il n'y a pas assez d'arguments dans la commande
  {
    printf("Veuillez préciser le chemin du fichier de matrix market.\nUtilisation: ./calcul_elements_propres chemin_fichier taille_sous_espace\n");
    return 1;
  }

  else // s'il y a assez d'arguments dans la commande
  {
    char *verifie_strtol; // pointeur utilisé par strtol pour vérifier que la taille spécifiée est bien un nombre
    size_t taille_sous_espace = strtol(argv[2], &verifie_strtol, 10);
    if (*verifie_strtol != '\0')
    {
      printf("La taille spécifiée pour le sous-espace n'est pas un nombre\nUtilisation: ./calcul_elements_propres chemin_fichier taille_sous_espace\n");
      return 1;
    }

    // initialisation des matrices et du vecteur yk
    gsl_spmatrix *A = lit_fichier_mat(argv[1]); // on importe la matrice du fichier
    gsl_matrix *B0 = gsl_matrix_calloc(taille_sous_espace, taille_sous_espace);
    gsl_matrix *B1 = gsl_matrix_calloc(taille_sous_espace, taille_sous_espace);
    gsl_matrix *Vm = gsl_matrix_calloc(A->size1, taille_sous_espace);
    gsl_vector *yk = gsl_vector_calloc(A->size1);

    gsl_vector_set_zero(yk); // le 1er élément de yk est égal à 1. Les autres sont égal à 0
    gsl_vector_set(yk, 0, 1);

    projection(A, B0, B1, Vm, yk, taille_sous_espace); // remplit B0, B1 et Vm et les affiche

    // TODO: yk=une combinaison linéaire des vecteurs propres au redémarrage

    // calcul et affichage de Em
    gsl_matrix *Em = inverse_matrix(B0);
    printf("EM:\n");
    affiche_matrice(Em);

    // calcul et affichage de Fm
    gsl_matrix *Fm = multiplie_matrices(Em, B1);
    printf("Fm:\n");
    affiche_matrice(Fm);

    // calcul des valeurs/vecteurs propres de Fm
    gsl_vector *valeurs_propres = gsl_vector_calloc(taille_sous_espace);                     // on stocke les valeurs propres dans un gsl_vector
    gsl_matrix *vecteurs_propres = gsl_matrix_alloc(taille_sous_espace, taille_sous_espace); // on stocke les vecteurs propres dans une gsl_matrix
    printf("matrice vecteurs_propres:\n");
    calcule_valeurs_et_vecteurs_propre(Fm, valeurs_propres, vecteurs_propres);

    printf("vecteurs_propres size1 = %ld et Vm size1 = %ld\n", vecteurs_propres->size1, Vm->size1);
    printf("vecteurs_propres size2 = %ld et Vm size2 = %ld\n", vecteurs_propres->size2, Vm->size2);

    // Partie 3
    // qi = Vm x ui avec ui = vecteur propre de Fm
    printf("\nProduit matrice-vecteur :\n");
    gsl_matrix *qi = gsl_matrix_alloc(A->size1, taille_sous_espace); // matrice qui contient les vecteurs qi
    gsl_vector *result = gsl_vector_alloc(A->size1);

    // calculs des vecteurs qi qu'on stocke dans la gsl_matrix qi
    for (size_t i = 0; i < taille_sous_espace; i++)
    {
      gsl_vector_view tempui = gsl_matrix_column(vecteurs_propres, i); // vecteurs propres d'une valeur propre (tt la colonne)
      gsl_vector *ui = &tempui.vector;
      // affiche_vecteur(ui, taille_sous_espace);
      produit_matrice_vecteur(Vm, ui, result);
      gsl_matrix_set_col(qi, i, result);

      /*for (size_t j = 0; j < vecteurs_propres->size2; j++)
      {
        gsl_vector *ui1 = gsl_vector_calloc(1);
        gsl_vector_set(ui1, 0, gsl_vector_get(ui, j));
        produit_matrice_vecteur(Vm, ui1);
      }*/
    }

    affiche_matrice(qi);

    // libération de la mémoire allouée
    gsl_spmatrix_free(A);
    gsl_matrix_free(B0);
    gsl_matrix_free(B1);
    gsl_matrix_free(Em);
    gsl_matrix_free(Fm);
    gsl_matrix_free(Vm);
    gsl_vector_free(yk);
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