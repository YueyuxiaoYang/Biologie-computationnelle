# Biologie-computationnelle
Projet Sam Mayer
https://github.com/bilalelhoudaigui/TCDS-v2

## Lancer une simulation
Dans le fichier test_simulation.py il faut choisir les paramètres suivants:

Le génome initial doit avoir le format gff et il faut changer le nom du fichier dans "input Genome information".

ratio= Pour espécifier le ration des insertions/deletions / inversions il faut. 0.1 signifie qu'il y a 1 insertion/deletion sur 10 inversions. 
T0= est la température pour l'algorithme de Montecarlo.
idl= est la longueur 


## mise à jour du 28 novembre 
## à faire 
1. save_tr_nbr does not change after modification of the Genome
1. sometimes simulation 'input' can not read all the genes.
1. visualize fitness to trace fitness value changing according to time
1. paralell computing for more than 1 genome at a time


####le 19 décembre 
#simulations à faire : 
  --> des résultats avec les paramètres initiaux : T=0.1, len=60 nucléotides,  ratio=1, nbre ARNpol=6
  --> simulation qui change le nbre d'ARNpol

#analyser les résultats : 
  --> taille des modifications 
  --> nbre ARNpol 

#choix des graphes (4 attention) : un graphe avec paramètres initiaux (effet modifications chromosomiques et fitness faible), le génome final ,mettre un graphe avec les paramètres optimaux, un graphe avec plusieurs graphes combinés.


#écire la définition de la fonction inversion 




#### le 11 janvier 
#simulations lancées : température et nombre ARNpoly : faire les graphes correspondants 
#écrire la définition de la fonction inversion 
#faire le graphe du génome final, graphe avec les paramètres optimaux 
