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
#include <omp.h>
#include <time.h>
#include <unistd.h>

// opérations sur vecteurs
double produit_scalaire(gsl_vector *yk, gsl_vector *yk_suivant);
double calcule_norme(gsl_vector *vecteur);
void produit_constante_vecteur(double constante, gsl_vector *vecteur, gsl_vector *resultat);
void soustrait_vecteur2_au_vecteur1(gsl_vector *vecteur1, gsl_vector *vecteur2, gsl_vector *resultat);

// opérations sur matrices
void produit_matrice_vecteur(gsl_matrix *m, gsl_vector *v, gsl_vector *resultat);
void produit_spmatrice_vecteur(gsl_spmatrix *m, gsl_vector *v, gsl_vector *resultat);
void calcule_valeurs_et_vecteurs_propre(gsl_matrix *matrix, gsl_vector *valeurs_propres, gsl_matrix *vecteurs_propres);
gsl_matrix *inverse_matrix(gsl_matrix *A);
gsl_matrix *multiplie_matrices(gsl_matrix *matrice1, gsl_matrix *matrice2);
void calcule_qi(gsl_spmatrix *A, gsl_matrix *qi, gsl_matrix *vecteurs_propres, gsl_matrix *Vm, size_t taille_sous_espace);
void projection(gsl_spmatrix *A, gsl_matrix *B0, gsl_matrix *B1, gsl_matrix *Vm, gsl_vector *yk, size_t taille_sous_espace);

// affichage
void affiche_matrice(gsl_matrix *A);
void affiche_vecteur(gsl_vector *v, int taille);

// importe une sp_matrix depuis un fichier du site matrix market
gsl_spmatrix *lit_fichier_mat(char nomFichier[]);

// renvoie 1 si l'epsilon est atteint en fin d'itération
int verifie_si_precision_atteinte(gsl_spmatrix *A, gsl_vector *valeurs_propres, gsl_matrix *qi, double epsilon);