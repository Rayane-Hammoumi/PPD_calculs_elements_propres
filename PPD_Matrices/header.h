#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <math.h>

void affiche_matrice(gsl_matrix *A);
gsl_spmatrix *lit_fichier_mat(char nomFichier[]);
void calcule_valeurs_propre(gsl_matrix *matrix);
gsl_matrix *inverse_matrix(gsl_matrix *A);
double produit_scalaire(gsl_vector *yk, gsl_vector *yk_suivant);
void projection(gsl_spmatrix *A, gsl_matrix *B, gsl_vector *yk);

#define taille_sous_espace_en_pourcentage 4 // par rapport à la matrice initiale A
// TODO:#define precision
