#include "header.h"

int main(int argc, char *argv[])
{
  if (argc == 1)
  {
    printf("Veuillez préciser le chemin du fichier de matrix market.\nutilisation: ./calcul_elements_propres cheminfichier\n");
    return 1;
  }

  else
  {
    // on importe la matrice du fichier
    gsl_spmatrix *A = lit_fichier_mat(argv[1]);

    // on définit le nombre de valeurs propres à calculer (taille du sous_espace)
    int taille_sous_espace_int = A->size1;
    taille_sous_espace_int = taille_sous_espace_int * taille_sous_espace_en_pourcentage / 100;
    size_t taille_sous_espace = taille_sous_espace_int;

    // initialisation de la matrice B et du vecteur yk
    gsl_matrix *B = gsl_matrix_calloc(taille_sous_espace, taille_sous_espace);
    gsl_vector *yk = gsl_vector_calloc(A->size1);

    gsl_vector_set_zero(yk); // le 1er élément de yk est égal à 1. Les autres sont égal à 0
    yk->data[0] = 1;

    projection(A, B, yk); // remplit la matrice B et l'affiche

    // TODO: yk=une combinaison linéaire des vecteurs propres au redémarrage

    // libération de la mémoire allouée
    gsl_spmatrix_free(A);
    gsl_matrix_free(B);
    gsl_vector_free(yk);
  }

  return 0;
}