scores_init <- function() {
  data.frame(
    date = c(),
    joueur_A1 = c(),
    joueur_A2 = c(),
    joueur_B1 = c(),
    joueur_B2 = c(),
    win = c(),
    game_length = c()
  )
}

add_scores <- function(daily_scores, scores, date) {
  #ordonner les 1vs 1 en premier pcq ca peut diminuer le temps de calcul
  #ne pas ordonner tous scores car on veut garder la temporalité entre les dates différentes.
  daily_scores <- rbind(daily_scores[is.na(daily_scores$joueur_A1), ], daily_scores[!is.na(daily_scores$joueur_A1), ])
  if(any(!apply(daily_scores[, paste0("joueur_", c("A1", "A2", "B1", "B2"))], 1, function(x) length(unique(x[!is.na(x)])) %in% c(2, 4)))) stop("Erreur dans le jeu de données")
  res <- rbind(scores, cbind(
    date = date,
    daily_scores
  ))
  dimnames(res)[[1]] <- 1:(nrow(scores) + nrow(daily_scores))
  res
}

