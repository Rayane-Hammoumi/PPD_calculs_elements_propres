#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

void affiche_matrice(gsl_matrix *A);
gsl_spmatrix *lit_fichier_mat(char nomFichier[]);
void calcule_valeurs_propre(gsl_matrix *matrix);
gsl_matrix *inverse_matrix(gsl_matrix *A);