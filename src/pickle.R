generate_likelihood_estimates_pickle <- function(dom_p = seq(0, 0.5, 0.01),
                                                 n_games = 100,
                                                 games_g = c(7, 11),
                                                 scores) {
  
  stopifnot(max(dom_p) <= 0.5)
  bootstrap <- function(simuls, n, dom) {
    replicate(n, {
      simuls_boot <- simuls[, sample(1:ncol(simuls), ncol(simuls), replace = TRUE)]
      apply(dom, 1, function(x) mean(simuls_boot[1, ] == x[1] & simuls_boot[2, ] == x[2]))
    })
  }
  
  res <- lapply(games_g, function(g) {
    dom_len_max <- max(scores[scores$game_len == g, c("score_A", "score_B")])
    possible_scores <- rbind(
      matrix(
        c(0:(g-2), rep(g, g-1), rep(g, g-1), c((g-2):0)),
        ncol = 2
      ),
      matrix(
        c(
          (g-1):(dom_len_max-1),
          (g+1):(dom_len_max+1),
          (g+1):(dom_len_max+1),
          (g-1):(dom_len_max-1)
        ),
        ncol = 2
      )
    )
    
    matrix(NA, nrow = length(dom_p), ncol = nrow(possible_scores),
           dimnames = list(dom_p, apply(possible_scores, 1, paste0, collapse = "-")))
  })
  
  for(i2 in seq_along(games_g)) {
    for(i1 in seq_along(dom_p)) {
      
      p <- dom_p[i1]
      g <- games_g[i2]
      dom_len_max <- max(scores[scores$game_len == g, c("score_A", "score_B")])
      
      simuls <- replicate(n_games, {
        scoreA <- 0
        scoreB <- 0
        serveA <- sample(c(TRUE, FALSE), 1)
        while(max(scoreA, scoreB) < g | abs(scoreA - scoreB) < 2) {
          win <- (1:2)[sample.int(2, 1, FALSE, c(p, 1-p))]
          if (win == 1) {
            if (serveA) scoreA <- scoreA + 1
            else serveA <- TRUE
          } else {
            if (serveA) serveA <- FALSE
            else scoreB <- scoreB + 1
          }
        }
        c(scoreA, scoreB)
      })
      
      possible_scores <- rbind(
        matrix(
          c(0:(g-2), rep(g, g-1), rep(g, g-1), c((g-2):0)),
          ncol = 2
        ),
        matrix(
          c(
            (g-1):(dom_len_max-1),
            (g+1):(dom_len_max+1),
            (g+1):(dom_len_max+1),
            (g-1):(dom_len_max-1)
          ),
          ncol = 2
        )
      )
      
      possible_scores <- cbind(
        possible_scores,
        apply(possible_scores, 1, function(x) mean(simuls[1, ] == x[1] & simuls[2, ] == x[2]))
      )
      
      #ggplot()+
      #  geom_point(aes(x = possible_scores[, 1], y = possible_scores[, 2], col = possible_scores[, 3]))
      
      
      # dist_boot <- bootstrap(simuls, 100, possible_scores[, 1:2])
      # dist_bound <- apply(dist_boot, 1, quantile, prob = c(0.1, 0.9))
      
      # dist <- sapply(0:(g-1), function(x) mean(simuls[2, ] == x))
      # ggplot()+
      #   geom_line(aes(x=0:(g-1), y=dist), col = "red")+
      #   geom_line(aes(x=0:(g-1), y=dnbinom(0:(g-1), g, p)))+
      #   geom_line(aes(x=0:(g-1), y=dist_bound[1, ]), linetype = "dashed")+
      #   geom_line(aes(x=0:(g-1), y=dist_bound[2, ]), linetype = "dashed")
      res[[i2]][i1, ] <- possible_scores[, 3]
    }
    print(paste0(i2, "/", length(games_g), collapse = ""))
  }
  names(res) <- games_g
  res
}

print("Loading...")
pickle_estim <- generate_likelihood_estimates_pickle(
  dom_p = seq(0, 0.5, 0.01),
  n_games = 3000,
  games_g = unique(as.numeric(scores$game_len)),
  scores
)

p_win_exact_not_vec_pickle <- function(
    p, scoreA, scoreB, game_len, win, pickle_estim
) {
  stopifnot((scoreA > scoreB & win) | (scoreA < scoreB & !win))
  
  dom_p <- as.numeric(rownames(pickle_estim[[as.character(game_len)]]))
  sapply(p, function(p) {
    # TODO pourrait arriver qu'on joue une partie sans écarts, dans ce cas,
    # il faudrait ajouter la donnée dans le generate avec un if
    if (p > max(dom_p)) {
      i <- which.min(abs(p - dom_p))
      pickle_estim[[as.character(game_len)]][
        i,
        paste0(scoreA, "-", scoreB, collapse = "")
      ]
    } else {
      i <- which.min(abs(1-p - dom_p))
      pickle_estim[[as.character(game_len)]][
        i,
        paste0(scoreB, "-", scoreA, collapse = "")
      ]
    }
  })
}

p_win_game_of_not_vec_pickle <- function(p, g, pickle_estim) {
  dom_p <- as.numeric(rownames(pickle_estim[[as.character(g)]]))
  tmp <- pickle_estim[[as.character(g)]]
  if (p > max(dom_p)) {
    1 - sum(tmp[which.min(abs(1-p - dom_p)), ])
  } else {
    sum(tmp[which.min(abs(p - dom_p)), ])
  }
}