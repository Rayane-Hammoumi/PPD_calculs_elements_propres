#include "header.h"

gsl_spmatrix *lit_fichier_mat(char nomFichier[])
{
    FILE *fp;
    gsl_spmatrix *A;

    fp = fopen(nomFichier, "r");
    A = gsl_spmatrix_fscanf(fp);
    fclose(fp);

    // Do something with the matrix, e.g. print it
    gsl_spmatrix_fprintf(stdout, A, "%.g");

    return A;
}