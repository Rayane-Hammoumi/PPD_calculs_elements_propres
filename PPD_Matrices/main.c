#include "header.h"

int main(int argc, char *argv[])
{
  if (argc < 3)
  {
    printf("Veuillez préciser le chemin du fichier de matrix market.\nUtilisation: ./calcul_elements_propres chemin_fichier taille_sous_espace\n");
    return 1;
  }

  else
  {
    char *verifie_strtol;
    size_t taille_sous_espace = strtol(argv[2], &verifie_strtol, 10);
    if (*verifie_strtol != '\0')
    {
      printf("La taille spécifiée pour le sous-espace n'est pas un nombre\nUtilisation: ./calcul_elements_propres chemin_fichier taille_sous_espace\n");
      return 1;
    }
    printf("%ld", taille_sous_espace);

    // on importe la matrice du fichier
    gsl_spmatrix *A = lit_fichier_mat(argv[1]);

    gsl_matrix *B0 = gsl_matrix_calloc(taille_sous_espace, taille_sous_espace);
    gsl_matrix *B1 = gsl_matrix_calloc(taille_sous_espace, taille_sous_espace);
    gsl_matrix *Vm = gsl_matrix_calloc(A->size1, taille_sous_espace);
    gsl_vector *yk = gsl_vector_calloc(A->size1);

    gsl_vector_set_zero(yk); // le 1er élément de yk est égal à 1. Les autres sont égal à 0
    yk->data[0] = 1;

    projection(A, B0, B1, Vm, yk, taille_sous_espace); // remplit B0 et B1 et les affiche

    // TODO: yk=une combinaison linéaire des vecteurs propres au redémarrage
    // affiche_matrice(multiplier_matrice(B1, B1));

    // libération de la mémoire allouée
    gsl_spmatrix_free(A);
    gsl_matrix_free(B1);
    gsl_vector_free(yk);
  }

  return 0;
}