source("src/update_scores.R")

faster_sample.int <- function(n, prob) {
  .Internal(sample(n, 1, FALSE, prob))
}

generate_likelihood_estimates_pickle <- function(dom_p = seq(0, 0.5, 0.005),
                                                 n_games = 100,
                                                 games_g = c(7, 11),
                                                 scores,
                                                 accelerator = FALSE) {
  
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
      
      acc_step <- 3 # à ajuster dans if(accelerator)
      acc_step_1 <- acc_step - 1
      prob_vector <- rep(dbinom(acc_step:0, acc_step, p) / choose(acc_step, acc_step:0), choose(acc_step, acc_step:0))
      stopifnot(all.equal(sum(prob_vector), 1))
      
      simuls <- replicate(n_games, {
        scoreA <- 0
        scoreB <- 0
        
        serveA <- sample(c(TRUE, FALSE), 1, prob = c(p, 1 - p))
        # supposons que le premier point sert à déterminer qui commence avec le service 
        
        if (accelerator) {
          
          while(max(scoreA, scoreB) < (g - acc_step_1)) {
            win_case <- (1:(2^acc_step))[faster_sample.int((2^acc_step), prob_vector)]
            if (win_case == 1) {
              #ppp
              if (serveA) scoreA <- scoreA + 1
              scoreA <- scoreA + 2
              serveA <- TRUE
            } else if (win_case == 2) {
              #ppq
              if (serveA) scoreA <- scoreA + 1
              scoreA <- scoreA + 1
              serveA <- FALSE
            } else if (win_case == 3) {
              #pqp
              if (serveA) {
                scoreA <- scoreA + 1
              } else {
                serveA <- TRUE
              }
            } else if (win_case == 4) {
              #qpp
              if (!serveA) scoreB <- scoreB + 1
              scoreA <- scoreA + 1
            } else if (win_case == 5) {
              #pqq
              if (serveA) {
                scoreA <- scoreA + 1
              } else {
                scoreB <- scoreB + 1
                serveA <- FALSE
              }
            } else if (win_case == 6) {
              #qpq
              if (!serveA) {
                scoreB <- scoreB + 1
              }
              serveA <- FALSE
            } else if (win_case == 7) {
              #qqp
              scoreB <- scoreB + 1
              if (!serveA) {
                scoreB <- scoreB + 1
              }
              serveA <- TRUE
            } else {
              #qqq
              if(!serveA) scoreB <- scoreB + 1
              scoreB <- scoreB + 2
              serveA <- FALSE
            }
          }
        }
        
        # code version la plus simple
        while(max(scoreA, scoreB) < g | abs(scoreA - scoreB) < 2) {
          win <- (1:2)[faster_sample.int(2, c(p, 1 - p))]
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
      
      if (p == 0.5) {
        possible_scores <- cbind(
          possible_scores,
          apply(
            possible_scores,
            1,
            function(x) mean(simuls[1, ] == x[1] & simuls[2, ] == x[2] | simuls[1, ] == x[2] & simuls[2, ] == x[1])
          )
        )
      } else {
        possible_scores <- cbind(
          possible_scores,
          apply(
            possible_scores,
            1,
            function(x) mean(simuls[1, ] == x[1] & simuls[2, ] == x[2])
          )
        )
      }
      
      possible_scores[, 3] <- possible_scores[, 3] / sum(possible_scores[, 3])
      
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

set.seed(2024L)
if ("pickle_estim_15000.RDS" %in% list.files()) {
  print("loading old pickle_estim_15000.RDS")
  pickle_estim <- readRDS("pickle_estim_15000.RDS")
  
  # TODO pas parfait car on pourrait avoir une partie qui dépasse le max précédent décart de 2 points
  if (!all(names(pickle_estim) %in% unique(as.numeric(scores[scores$serve_for_pt, "game_len"])))) {
    new_g <- names(pickle_estim)[names(pickle_estim) %in% unique(as.numeric(scores[scores$serve_for_pt, "game_len"]))]
    pickle_estim2 <- generate_likelihood_estimates_pickle(
      dom_p = seq(0, 0.5, 0.005),
      n_games = 15000,
      games_g = new_g,
      scores
    )
    saveRDS(c(pickle_estim, pickle_estim2), "pickle_estim_15000.RDS")
  }
} else {
  print("Loading...")
  pickle_estim <- generate_likelihood_estimates_pickle(
    dom_p = seq(0, 0.5, 0.005),
    n_games = 15000,
    games_g = unique(as.numeric(scores[scores$serve_for_pt, "game_len"])),
    scores
  )
  saveRDS(pickle_estim, "pickle_estim_15000.RDS")
}

p_win_exact_not_vec_pickle <- function(
    p, scoreA, scoreB, game_len, win, pickle_estim
) {
  stopifnot((scoreA > scoreB & win) | (scoreA < scoreB & !win))
  if (win == 0) {
    p <- 1 - p
    tmp <- scoreB
    scoreB <- scoreA
    scoreA <- tmp
  }
  dom_p <- as.numeric(rownames(pickle_estim[[as.character(game_len)]]))
  sapply(p, function(p) {
    # TODO pourrait arriver qu'on joue une partie sans écarts, dans ce cas,
    # il faudrait ajouter la donnée dans le generate avec un if
    if (p <= max(dom_p)) {
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
  
  col_where_A_wins <- c(g:(g + g - 2), (ncol(tmp) - ((ncol(tmp) - (g-1) * 2) / 2)):ncol(tmp))
  
  if (p > max(dom_p)) {
    1 - sum(tmp[which.min(abs(1-p - dom_p)), col_where_A_wins])
  } else {
    sum(tmp[which.min(abs(p - dom_p)), col_where_A_wins])
  }
}

stopifnot(all.equal(
  p_win_exact_not_vec_pickle(0.4, 11, 7, 11, 1, pickle_estim),
  p_win_exact_not_vec_pickle(1 - 0.4, 7, 11, 11, 0, pickle_estim)
))
stopifnot(all.equal(
  p_win_exact_not_vec_pickle(0.7, 11, 7, 11, 1, pickle_estim),
  p_win_exact_not_vec_pickle(1 - 0.7, 7, 11, 11, 0, pickle_estim)
))
stopifnot(all.equal(
  p_win_exact_not_vec_pickle(0.5, 11, 7, 11, 1, pickle_estim),
  p_win_exact_not_vec_pickle(0.5, 7, 11, 11, 0, pickle_estim)
))

stopifnot(all.equal(
  likelihood_1vs1_exact(joint_density,
                        likelihood_1vs1_exact_prob_win_1_pt(joint_density, c("A", "B")),
                        score = c(date = NA, joueur_A1 = NA, joueur_A2 = "A", joueur_B1 = "B", joueur_B2 = NA, win = 1, score_A = 11, score_B = 7, game_len = 11, serve_for_pt = TRUE),
                        dataset = "pickle"),
  likelihood_1vs1_exact(joint_density,
                        likelihood_1vs1_exact_prob_win_1_pt(joint_density, c("B", "A")),
                        score = c(date = NA, joueur_A1 = NA, joueur_A2 = "B", joueur_B1 = "A", joueur_B2 = NA, win = 0, score_A = 7, score_B = 11, game_len = 11, serve_for_pt = TRUE),
                        dataset = "pickle")
))

