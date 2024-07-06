Pour mettre à jour le csv : aller dans Excel, télécharger une copie du xlsx dans ordi. Ensuite, l'ouvrir et enregistrer sous format CSV UTF-8

# TODO
- Améliorer la considération des scores exacts si non crédible. Surement pas mettre un seuil fixe comme actuellement, y aller progressivement. 
- Il faudrait aussi s'assurer que quelquun de bas, ne puisse pas arrêter de jouer pour remonter à 50. 
- Ajouter une manière de pénaliser les gens vraiment pas crédible? Ou modifier le calcul de prob dans le ELO en tenant compte de la variance
- maybe calibrer la priori selon la distribution empirique du monde après plusieurs mois. Possiblement mettre une masse à gauche pour résoudre le problème de drift (incitatif aux faibles à ne plus jouer pour monter)
- mettre tous les hp dans un fichier de config R, incluant les distributions priori
- faire une graphe avec des arêtes pour montré qui a joué contre qui.
- Possiblement faire un graphique à 2 axe =, coin haut droit est 100, coin bas gauche est 0, et il y aurait des cercles.
- Modifier le drift, de manières à ajouter des probs juste aux points dans le domaine et voisins au domaine p>0 du joueur, avec les probs sommées, pour pas trop ogmenter la dimensionnalité
- L'ajout de plusieurs joueurs tard dans le processus augmente bcp la dimensionnalité et parfois
dépasse le absolute_max_dim donc on perd beaucoup de probs.
- Essayer de modéliser les distributions par des splines de degré 1-3 mettons. Chq personne aurait peu de params.