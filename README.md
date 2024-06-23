Pour mettre à jour le csv : aller dans Excel, télécharger une copie du xlsx dans ordi. Ensuite, l'ouvrir et enregistrer sous format CSV UTF-8

# TODO
- ajouter une matrice qui montre le nb de parties jouées entre tout le monde
- Ajouter un indice de crédibilité selon la distribution à posteriori
- Ajouter une fonction de recommandation des prochaines parties qu'il faut obtenir pour stabiliser les scores.
- débugger le bug quand include_exact_points <- T
- améliorer le kernel smoothing en ajustant la largeur des kernels et retirant le capping
- Bug Jason premier
- Pour les personnes non crédibles, utiliser le score exact par défaut?
- changer le include_exact_points <- T pour juste changer la formule de likelihood pour tenir comte du score exact. ca va etre plus rapide et plus précis pour éviter que l'ordonnancement des parties
ait un impact.