#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <ctype.h>

void affiche_matrice(gsl_matrix *A);
void affiche_vecteur(gsl_vector *v, int taille);
gsl_spmatrix *lit_fichier_mat(char nomFichier[]);
void calcule_valeurs_et_vecteurs_propre(gsl_matrix *matrix, gsl_vector* valeurs_propres, gsl_matrix * vecteurs_propres);
gsl_matrix *inverse_matrix(gsl_matrix *A);
double produit_scalaire(gsl_vector *yk, gsl_vector *yk_suivant);
void projection(gsl_spmatrix *A, gsl_matrix *B0, gsl_matrix *B1, gsl_matrix *Vm, gsl_vector *yk, size_t taille_sous_espace);
gsl_matrix *multiplie_matrices(gsl_matrix *matrice1, gsl_matrix *matrice2);
void remplis_colonne_matrice_avec_vecteur(gsl_matrix *matrice, gsl_vector *vecteur, size_t numero_colonne);
gsl_vector* produit_matrice_vecteur(gsl_matrix *m, gsl_vector *v);

// TODO:#define precision
