#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

gsl_matrix * lit_fichier_mat(char nomFichier[])
{
    int m, n;
    gsl_matrix *A;
    FILE *f;
    char line[1024];

    f = fopen(nomFichier, "r");
    if (f == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }
    // on ignore les lignes ou il y a du texte
    while (fgets(line, sizeof(line), f)) {
        if (sscanf(line, "%d %d", &m, &n) == 2) {
            //Found the matrix size
            break;
        }
    }
    printf("%d %d\n",m, n);
    A = gsl_matrix_alloc(m, n);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double value;
            fscanf(f, "%lf", &value);
            gsl_matrix_set(A, i, j, value);
        }
    }
    fclose(f);
    //Afficher la matrice
    gsl_matrix_fprintf (stdout, A, "%g");
    return A;
}

int main(int argc, char*argv[])
{
  /*int n = 3;
  gsl_matrix *A = gsl_matrix_alloc(n, n);
  // on remplit la matrice A avec des valeurs
  
  gsl_matrix_set(A, 0, 0, 2.0);
  gsl_matrix_set(A, 0, 1, 1.0);
  gsl_matrix_set(A, 0, 2, 0.0);
  gsl_matrix_set(A, 1, 0, 1.0);
  gsl_matrix_set(A, 1, 1, 2.0);
  gsl_matrix_set(A, 1, 2, 1.0);
  gsl_matrix_set(A, 2, 0, 0.0);
  gsl_matrix_set(A, 2, 1, 1.0);
  gsl_matrix_set(A, 2, 2, 2.0);*/
  
  if (argc > 1)
  {
    gsl_matrix *B = lit_fichier_mat(argv[1]);
    int n = B->size1;

    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
    gsl_vector *eval = gsl_vector_alloc(n);
    gsl_matrix *evec = gsl_matrix_alloc(n, n);
    //gsl_eigen_symmv(A, eval, evec, w);
    gsl_eigen_symmv(B, eval, evec, w);
    gsl_eigen_symmv_free(w);

    // afficher les valeurs propres 
    for (int i = 0; i < n; i++)
    {
      printf("valeur propre = %g\n", gsl_vector_get(eval, i));
    }

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    //gsl_matrix_free(A);
    gsl_matrix_free(B);
  }

  return 0;
}

