# Biologie-computationnelle
Projet Sam Mayer
https://github.com/bilalelhoudaigui/TCDS-v2

Le code simule l'évolution d'un génome en utilisant des insertions/délétions et inversions. Le but du code est de génerer un génome qui soit le mieux adapté à son environement au niveau de transcription des gènes. Le génome est modifié au cours du temps et la transcription des gènes est calculée en utilisant le package TCDS-v2. La transcription cible des gènes se trouve dans le fichier input1/environment.data. Il faut trouver une conformation des gènes qui ait une transcription proche de la cible. Avec ce code il est possible de tester plusieurs paramètres pour trouver la bonne conformation du génome. 

# Lancer une simulation

## Paramètres
Dans le fichier test_simulation.py il faut choisir les paramètres suivants:
  -->Le génome initial doit avoir le format gff et il faut changer le nom du fichier dans "input Genome information".

  --> ratio = Pour spécifier le ratio des (insertions/deletions)/inversions. Ce ratio doit être égal à 0.1, signifiant qu'il y a          1 insertion/deletion sur 10 inversions. 
  
  --> T0= est la température pour l'algorithme de Montecarlo.
  --> idl= est la longueur des insertions/deletions
  --> iter_num= est le nombre de générations du génome
  --> rep_num= est le nombre de répétitions d'un simulation

  --> INI_file= le nom du fichier params.ini du package TCDS.

Si vous voudriez tester plusieurs paramètres vous pouvez changer la variable T0_list avec les différentes valeurs à tester.
Dans ce cas là nous avons testé différentes valeurs de tempétature.

## Remarques 
Pour lancer le code il faut lancer le fichier test_simulations.py.

Les fichiers input gff et params.ini doivent être placés dans le même dossier du test_simulations.py.

## Fichiers du sortie

fit_total.output : fichier avec la fitness au cours des générations du génome.
mutation_total.output : fichier avec l'information des mutations et des inversions au cours des générations.
name_total.output : fichier avec l'information des paramètres de la simulation.

Pour modifier le nom des fichiers ou le nom du dossier de sortie vous pouvez les changer à la fin du code test_simulations.py.
