Afin de compiler ce projet, en plus d'OpenMP, vous aurez besoin de la bibliothèque GSL.

#Installation sous linux: sudo apt-get install libgsl-dev

#Documentation: https://www.gnu.org/software/gsl/doc/html/vectors.html#

#Compiler le projet: make

#Exécuter le projet: ./calcul_elements_propres chemin_fichier_sparse_matrice taille_sous_espace précision

NB: taille_sous_espace est un entier, précision est un nombre décimal (0.00000000001 par exemple)

Vous disposez de plusieurs fichiers .mtx dans ce dossier sur lequel vous pouvez exécuter le programme.