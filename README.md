Pour mettre à jour le csv : aller dans Excel, télécharger une copie du xlsx dans ordi. Ensuite, l'ouvrir et enregistrer sous format CSV UTF-8

# TODO
- Améliorer la considération des scores exacts si non crédible. Surement pas mettre un seuil fixe comme actuellement, y aller progressivement. 
- Il faudrait aussi s'assurer que quelquun de bas, ne puisse pas arrêter de jouer pour remonter à 50. 
- Ajouter une manière de pénaliser les gens vraiment pas crédible? Ou modifier le calcul de prob dans le ELO en tenant compte de la variance
- maybe calibrer la priori selon la distribution empirique du monde après plusieurs mois. Possiblement mettre une masse à gauche pour résoudre le problème de drift (incitatif aux faibles à ne plus jouer pour monter)
- mettre tous les hp dans un fichier de config R, incluant les distributions priori
- faire une graphe avec des arêtes pour montré qui a joué contre qui.
- rescale the skill parameter between 0 and 1
- ajouter densité avec écart de 2 points ou non 
- réduire le hp de drift
- pour l'initiatisation des theta dans l'optimisation, utiliser le résultat précédent en mémoire pour accélérer la convergence.
- Possiblement faire un graphique à 2 axe =, coin haut droit est 100, coin bas gauche est 0, et il y aurait des cercles.
- Dans le graphique du ranking dans le temps, juste mettre des points quand les gens ont joué des parties, lignes si c'est juste du drift ou des ajustements par rapport aux autres.
- Ajouter la crédibilité dans le graph de ranking actuel, pt juste avec la hauteur des barres.