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

scores_per_pt_converter <- function(score) {
  A <- as.numeric(score["score_A"])
  B <- as.numeric(score["score_B"])
  res <- as.data.frame(matrix(score, nrow = A + B, ncol = 8, byrow=T, dimnames = list(NULL, names(score))))
  if(A >= 1) {
    res[1:A, "win"] <- 1
    res[1:A, "score_A"] <- 1
    res[1:A, "score_B"] <- 0
  }
  if(B >= 1) {
    res[(A + 1) : (A + B), "win"] <- 0
    res[(A + 1) : (A + B), "score_A"] <- 0
    res[(A + 1) : (A + B), "score_B"] <- 1
  }
  res$game_length <- 1
  res
}

scores_stats <- function(scores, players) {
  #games played per players
  games_played <- sapply(names(players), function(nom) sum(apply(scores[, paste0("joueur_", c("A1", "A2", "B1", "B2"))], 1, function(game) nom %in% game)))
  
  list(
    games_played = games_played[order(games_played, decreasing = T)],
    
    #points gained/loss per players
    pts = sapply(names(sort(games_played, decreasing = T)), function(nom) {
      keep_A <- apply(scores[, paste0("joueur_", c("A1", "A2"))], 1, function(game) nom %in% game)
      keep_B <- apply(scores[, paste0("joueur_", c("B1", "B2"))], 1, function(game) nom %in% game)
      pts_win <- sum(as.numeric(scores[keep_A, "score_A"])) + sum(as.numeric(scores[keep_B, "score_B"]))
      pts_played <- sum(as.numeric(scores[keep_A | keep_B, "score_A"])) + sum(as.numeric(scores[keep_A | keep_B, "score_B"]))
      
      c(pts_win = round(pts_win, 1), pts_played = round(pts_played, 1), ratio = round(pts_win / pts_played, 3))
    }
    )
  )
}

games_matchups <- function(scores, players) {
  noms <- names(players)
  combins <- combn(noms, 2)
  games_1vs1 <- apply(scores, 1, function(x) is.na(x["joueur_A1"]))
  res <- apply(combins, 2, function(noms) sum(apply(scores[games_1vs1, c("joueur_A2", "joueur_B1"), drop=F], 1, function(x) all(noms %in% x))))
  data.frame(A=combins[1, ][order(res, decreasing = T)], B=combins[2, ][order(res, decreasing = T)], count=res[order(res, decreasing = T)])
}
