recommend_next_game <- function(players, names_present = NULL) {
  
  if(is.null(names_present)) {
    noms <- names(players)
  } else {
    noms <- names_present
  }
  combins <- combn(noms, 2)
  
  info_gained <- apply(combins, 2, function(pair) {
    new_dists <- update_scores(
      players[pair],
      #fake 2 games, one won, other lost, to see how much it would change distributions
      scores=data.frame(
        date=as.character(Sys.Date()),
        "joueur_A1"=NA, "joueur_A2"=pair[1], "joueur_B1"=pair[2], "joueur_B2"=NA,
        win=c(1, 0), score_A=c(5, 3), score_B=c(3, 5), game_len=5
      )
    )
    
    mean(mapply(wassertein, new_dists, players[pair]))
  })
  
  rbind(combins, round(info_gained, 2))[, head(order(info_gained, decreasing = T), 10)]
}
