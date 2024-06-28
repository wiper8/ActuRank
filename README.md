Pour mettre à jour le csv : aller dans Excel, télécharger une copie du xlsx dans ordi. Ensuite, l'ouvrir et enregistrer sous format CSV UTF-8

# TODO
- ajouter une matrice qui montre le nb de parties jouées entre tout le monde
- weights dans Likelihood : retrouver une manière de weigther les games récentes
- Améliorer la considération des scores exacts si non crédible. Surement pas mettre un seuil fixe comme actuellement, y aller progressivement. Peut-être ajouter le drift à tout le monde. Il faudrait aussi s'assurer que quelquun de bas, ne puisse pas arrêter de jouer pour remonter à 50. 
- Ajouter une manière de pénaliser les gens vraiment pas crédible? Ou modifier le calcul de prob dans le ELO en tenant compte de la variance
- maybe calibrer la priori selon la distribution empirique du monde après plusieurs mois. Possiblement mettre une masse à gauche pour résoudre le problème de drift (incitatif aux faibles à ne plus jouer pour monter)
- mettre tous les hp dans un fichier de config R, incluant les distributions priori
- faire une graphe avec des arêtes pour montré qui a joué contre qui.