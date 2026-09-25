source("src/update_scores.R")

recommend_next_game <- function(players, names_present = NULL) {
  
  if(is.null(names_present)) {
    noms <- names(players)
  } else {
    noms <- names_present
  }
  combins <- combn(noms, 2)
  
  info_gained <- apply(combins, 2, function(pair) {
    
    ranks <- sapply(players[pair], function(distr) calculate_skill(distr, players))
    ranks <- sort(ranks, decreasing = T)
    
    new_dists <- update_scores_exact(
      players[pair],
      # fake game, one won, to see how much it would change distributions if the current 
      # best wins it
      scores=data.frame(
        date=as.character(Sys.Date()),
        "joueur_A1"=NA, "joueur_A2"=pair[1], "joueur_B1"=pair[2], "joueur_B2"=NA,
        win=1, score_A=5, score_B=3, game_len=5
      )
    )
    
    mean(mapply(wassertein, new_dists, players[pair]))
  })
  
  rbind(combins, round(info_gained, 2))[, head(order(info_gained, decreasing = T), 10)]
}

recommend_fair_teams <- function(joint_density, present_players = NULL) {
  if (!is.null(present_players)) {
    # filtrer pour les joueurs présents
    joint_density$names <- present_players
    joint_density$domains <- joint_density$domains[present_players]
    joint_density$joint_distr <- joint_density$joint_distr[, c(present_players, "p")]
    joint_density$grid_id <- joint_density$grid_id[, present_players]
  }

  pairs <- as.matrix(pair_vs_pair(length(joint_density$names)))
  combins <- joint_density$names[pairs] |> matrix(nrow = nrow(pairs))
  probs_best_team <- apply(combins, 1, function(comb) {
    sum((abs(likelihood_2vs2_exact_prob_win_1_pt(joint_density, comb) - 0.5) + 0.5) *
      joint_density$joint_distr$p)
  })
  cbind(probs_best_team, combins)[order(probs_best_team), ]
  res <- combins[which.min(probs_best_team), ]
  paste0(res[1], " & ", res[2], " vs ", res[3], " & ", res[4])
}

# code chatGPT
pair_vs_pair <- function(n) {
  players <- seq_len(n)
  
  groups <- combn(players, 4, simplify = FALSE)
  
  do.call(rbind, lapply(groups, function(x) {
    rbind(
      c(x[1], x[2], x[3], x[4]),
      c(x[1], x[3], x[2], x[4]),
      c(x[1], x[4], x[2], x[3])
    )
  })) |>
    as.data.frame() |>
    setNames(c("team1_p1", "team1_p2", "team2_p1", "team2_p2"))
}
