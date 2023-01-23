#include "header.h"

gsl_matrix *lit_fichier_mat(char nomFichier[])
{
    int m, n;
    gsl_matrix *A;
    FILE *f;
    char line[1024];

    f = fopen(nomFichier, "r");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    // on ignore les lignes ou il y a du texte
    while (fgets(line, sizeof(line), f))
    {
        if (sscanf(line, "%d %d", &m, &n) == 2)
        {
            // Found the matrix size
            break;
        }
    }
    printf("%d %d\n", m, n);
    A = gsl_matrix_alloc(m, n);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double value;
            fscanf(f, "%lf", &value);
            gsl_matrix_set(A, i, j, value);
        }
    }
    fclose(f);
    // Afficher la matrice
    gsl_matrix_fprintf(stdout, A, "%g");
    return A;
}