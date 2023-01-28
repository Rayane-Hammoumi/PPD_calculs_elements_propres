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
    yk->data[0] = 1;

    projection(A, B0, B1, Vm, yk, taille_sous_espace); // remplit B0, B1 et Vm et les affiche

    // TODO: yk=une combinaison linéaire des vecteurs propres au redémarrage

    // calcul et affichage de Em
    gsl_matrix *Em = inverse_matrix(B0);
    printf("EM:\n");
    affiche_matrice(Em);

    // calcul et affichage de Fm
    // gsl_matrix *Fm = multiplie_matrices(Em, B1);
    gsl_matrix *Fm = multiplie_matrices(Em, B1);
    printf("Fm:\n");
    affiche_matrice(Fm);

    // libération de la mémoire allouée
    gsl_spmatrix_free(A);
    gsl_matrix_free(B0);
    gsl_matrix_free(B1);
    gsl_matrix_free(Em);
    gsl_matrix_free(Fm);
    gsl_vector_free(yk);
  }
  gsl_matrix *m = gsl_matrix_alloc(2, 2);
  gsl_matrix_set(m, 0, 0, 1.0);
  gsl_matrix_set(m, 0, 1, 2.0);
  gsl_matrix_set(m, 1, 0, 4.0);
  gsl_matrix_set(m, 1, 1, 3.0);
  
  affiche_matrice(m);

  gsl_vector* valeurs_propres = gsl_vector_calloc(2);
  gsl_vector_set_zero(valeurs_propres);
  gsl_matrix * vecteurs_propres = gsl_matrix_alloc(2, 2);

  
  calcule_valeurs_propre(m, valeurs_propres, vecteurs_propres);

  return 0;
}