library(ggplot2)
library(igraph)

show_skill_level <- function(players) {
  ranks <- sapply(players, function(distr) calculate_skill(distr, players))
  ranks <- ranks[order(ranks, decreasing = TRUE)]
  
  x <- seq(0, 1, 0.01)
  
  graph_data <- data.frame(
    skill = ranks,
    player = names(ranks)
  )
  
  print(ggplot()+
          geom_vline(aes(xintercept = skill, col = player), data=graph_data, linewidth=2)+
          scale_color_discrete(breaks = graph_data$player[order(graph_data$skill, decreasing = T)])+
          theme_bw()+xlab("Skill")+
          xlim(0, 1)
  )
  ranks
}

smooth_distr <- function(distr, step = 5) {
  stopifnot(100 %% step == 0)
  res <- cbind(mu = seq(0 + step / 2, 100 - step / 2, step), 
               p = sapply(seq(0 + step / 2, 100 - step / 2, step),
                          function(x_seuil) sum(distr[, "p"][distr[, "mu"] <= (x_seuil + step / 2) & distr[, "mu"] > (x_seuil - step / 2)])))
  res <- res[res[, 2] > 0, , drop = FALSE]
  if(abs(sum(res[, "p"]) - 1) > 0.02) stop("bug de somme de prob trop loin de 1")
  res[, "p"] <- res[, "p"] / sum(res[, "p"])
  res
}

complement_credibilite <- function(players) {
  distr <- do.call(rbind, players)
  distr <- rbind(c(distr[1, 1], p=0), distr)
  x <- distr[, "mu"]
  y <- distr[, "p"]
  ordre <- order(x)
  y <- y[ordre]
  y <- y/sum(y)
  x <- x[ordre]
  
  smooth_distr(matrix(c(x, y), ncol=2, dimnames = list(NULL, c("mu", "p"))), step = 5)
}

credibilise <- function(distr, players, seuil = 0.7) {
  z <- min(1, compute_credibility(distr) / seuil)
  complement <- complement_credibilite(players)
  rbind(cbind(mu = distr[, "mu"], p = z * distr[, "p"]), cbind(mu = complement[, "mu"], p = (1 - z) * complement[, "p"]))
}

#TODO ajouter une version exacte?
show_current_probs <- function(players) {
  ranks <- sapply(players, function(distr) calculate_skill(distr, players))
  ordre <- order(ranks, decreasing = TRUE)
  players <- players[ordre]
  ranks <- ranks[ordre]
  pairs <- t(combn(1:length(ranks), 2))
  #players_credibilise <- lapply(players, function(distr) credibilise(distr, players))
  
  players <- lapply(players, simplifier_domain)
  
  prob <- mapply(function(distr_S1, distr_S2) {
    
    distr_S1_S2 <- distr_F1_F2_1vs1(distr_S1, distr_S2)
    
    MS1 <- rep(lapply(distr_S1[, "mu"] / 100, transition_matrix), nrow(distr_S2))
    MS2 <- rep(lapply(distr_S2[, "mu"] / 100, transition_matrix), each = nrow(distr_S1))
    
    
    sum(mapply(function(MS1, MS2) prob_win_point_1vs1_knowing_skills(MS1, MS2),
               MS1, MS2) * distr_S1_S2[, "p1"] * distr_S1_S2[, "p2"])
  },
  players[pairs[, 1]],#players_credibilise[pairs[, 1]],
  players[pairs[, 2]]#players_credibilise[pairs[, 2]]
  )
  
  data.frame(A=names(ranks)[pairs[, 1]], B=names(ranks)[pairs[, 2]],
             prob = round(prob, 3), prob_win_11 = round(sapply(prob, function(p) p_win_game_of(p, 11)), 3))
}

show_current_ranking <- function(players, scores, init_theta = NULL, show_credibility = FALSE) {
  probs <- show_current_probs(players)
  
  ranks <- sapply(players, function(distr) calculate_skill(distr, players))
  ordre <- order(ranks, decreasing = TRUE)
  players <- players[ordre]
  ranks <- ranks[ordre]
  init_theta <- init_theta[names(players)]
  pairs <- t(combn(1:length(ranks), 2))
  
  # Trouver les composantes fortement connexes d'un graphe
  pairings <- players_pairs(scores)
  
  pairings <- pairings[sapply(pairings, function(x) all(x %in% names(players)))]
  
  edges <- matrix(unlist(pairings), nrow=2)
  g <- graph(edges, directed = FALSE)
  clust <- components(g)$membership
  clust <- clust[names(players)]
  
  contraintes <- matrix(0, nrow = max(clust) * 2, ncol = length(players))
  
  for(i in 1:max(clust)) {
    contraintes[(i - 1) * 2 + 1, clust == i] <- 1
    contraintes[(i - 1) * 2 + 2, clust == i] <- -1
  }
  
  if (is.null(init_theta) | any(is.na(names(init_theta)))) {
    init_theta <- rep(50, length(players))
    for(i in 1:max(clust)) {
      init_theta[clust == i] <- qnorm(seq(0.1, 0.9, length.out = sum(clust == i)), 50, 15)
    }
  }
  
  
  #weighter un peu par la crédibilité
  weights <- sapply(players[names(ranks)], compute_credibility)
  weights <- mapply(sum, weights[pairs[, 1]], weights[pairs[, 2]])
  to_optim <- function(Forces) {
    estim <- 1 / (1 + 10^(-(Forces[pairs[, 1]] - Forces[pairs[, 2]]) / 20))
    mean(weights * (estim - probs[, "prob_win_11"])^2)
  }
  
  #min 0
  # res <- constrOptim(rep(50, length(players)), to_optim, grad = NULL,
  #             ui = rbind(diag(length(players))),
  #             ci = rep(0, length(players)))

  ui <- rbind(diag(length(players)), -diag(length(players)), rep(1, length(players)), rep(-1, length(players)), contraintes)
  ci <- c(rep(0, length(players)), rep(-100, length(players)), 49.5*length(players), -50.5*length(players), sapply(1:max(clust), function(i) c(49.5, -50.5) * sum(clust == i)))
  
  #enlever les erreurs de contrainte == 0
  if(any(ui %*% init_theta - ci <= 0)) {
    which_constraint_violated <- which(ui %*% init_theta - ci <= 0)
    for(i in which_constraint_violated) {
      init_theta[ui[i, ] != 0] <- pmin(pmax(0.1, init_theta[ui[i, ] != 0] + (50 - mean(init_theta[ui[i, ] != 0]))), 99.9)
    }
  }
  if(any(ui %*% init_theta - ci <= 0)) {
    init_theta <- rep(50, length(players))
  }
  
  #moyenne 50 et min 0 et max 100
  #pourrait donner des bugs si le monde est très dispersé
  res <- constrOptim(init_theta, to_optim, grad = NULL, ui = ui, ci = ci)

  
  #scores entre 1 et 100
  # res <- constrOptim(rep(50, length(players)), to_optim, grad = NULL,
  #                   ui = rbind(diag(length(players)), -diag(length(players))),
  #                   ci = c(rep(1, length(players)), rep(-100, length(players))))
  
  #sum(scores) = 50*n
  # res <- constrOptim(rep(50, length(players)), to_optim, grad = NULL,
  #             ui = rbind(diag(length(players)), -diag(length(players)), rep(1, length(players)), rep(-1, length(players))),
  #             ci = c(rep(1, length(players)), rep(-100, length(players)), 49.5*length(players), -50.5*length(players)))
  names(res$par) <- names(players)
  
  if(sqrt(res$value) > 0.02) warning("Convergence non parfaite des scores")
  
  credibl <- sapply(players[names(res$par)], compute_credibility)
  
  graph_data <- data.frame(
    score = res$par,
    player = names(res$par),
    credibility = credibl
  )
  
  if(show_credibility) {
    
    print(
      ggplot()+
        geom_linerange(aes(ymin = 0, ymax = credibility, x = score, col = player), data=graph_data, linewidth=2)+
        geom_text(aes(x=score, y = runif(nrow(graph_data), 0, credibility), label=player), data=graph_data, angle=45)+
        scale_color_discrete(breaks = graph_data$player[order(graph_data$score, decreasing = T)])+
        #xlim(0, 100)+
        theme_bw()+
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+ylab("")
    )
  } else {
    
    print(
      ggplot()+
        geom_vline(aes(xintercept = score, col = player), data=graph_data, linewidth=2)+
        geom_text(aes(x=score, y = sample(nrow(graph_data), nrow(graph_data)), label=player), data=graph_data, angle=45)+
        scale_color_discrete(breaks = graph_data$player[order(graph_data$score, decreasing = T)])+
        #xlim(0, 100)+
        theme_bw()+
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+ylab("")
    )
  }
  
  
  round(sort(res$par, decreasing = T), 1)
}


show_detailed_skill <- function(players) {
  
  ranks <- sapply(players, function(distr) calculate_skill(distr, players))
  ranks <- ranks[order(ranks, decreasing = TRUE)]
  print(ranks)
  
  graph_data <- data.frame(
    skill = ranks,
    player = names(ranks)
  )
  
  players2 <- lapply(players, function(distr) {
    distr[distr[, "p"] > 0, , drop = FALSE]
  })
  
  graph_data2 <- do.call(rbind, mapply(function(x, n) data.frame(x, "player"=n), players2, names(players2), SIMPLIFY = FALSE))
  
  tmp <- table(graph_data2$player)
  graph_data2_one_point <- graph_data2[graph_data2$player %in% names(tmp)[tmp == 1], ]
  
  ggplot()+
    #geom_vline(aes(xintercept = skill, col = player), data=graph_data, linewidth=2)+
    geom_segment(aes(x=mu/100, xend=mu/100, y = 0, yend = p, col = player), data = graph_data2_one_point,
               linewidth=1, alpha=0.8)+
    geom_line(aes(x=mu/100, y=p, col = player), data = graph_data2[!graph_data2$player %in% graph_data2_one_point$player, ],
              linewidth=1, alpha=0.8)+
    xlab("Skill")+ylab("Likelihood")+
    theme_bw()+
    xlim(0, 1)
  
}

show_detailed_skill_per_player <- function(players) {
  
  ranks <- sapply(players, function(distr) calculate_skill(distr, players))
  ranks <- ranks[order(ranks, decreasing = TRUE)]
  print(ranks)
  
  graph_data <- data.frame(
    skill = ranks,
    player = names(ranks)
  )
  
  players2 <- lapply(players, function(distr) {
    distr[distr[, "p"] > 0, , drop = FALSE]
  })
  
  graph_data2 <- do.call(
    rbind,
    mapply(
      function(x, n) data.frame(x, "player"=n), players2, names(players2),
      SIMPLIFY = FALSE
    )[names(sort(ranks, decreasing = TRUE))]
  )
  
  graph_data2$player <- factor(graph_data2$player, levels = unique(graph_data2$player))
  
  tmp <- table(graph_data2$player)
  graph_data2_one_point <- graph_data2[graph_data2$player %in% names(tmp)[tmp == 1], ]
  
  ggplot()+
    #geom_vline(aes(xintercept = skill, col = player), data=graph_data, linewidth=2)+
    geom_segment(aes(x=mu/100, xend=mu/100, y = 0, yend = p, col = player), data = graph_data2_one_point,
                 linewidth=1, alpha=0.8)+
    geom_line(aes(x=mu/100, y=p), data = graph_data2[!graph_data2$player %in% graph_data2_one_point$player, ],
              linewidth=1, alpha=0.8)+
    xlab("Skill")+ylab("Likelihood")+
    theme_bw()+
    xlim(0, 1)+
    facet_grid(rows = vars(player), scales = "free_y")+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "none")
}

show_IC_skill <- function(players) {
  
  ranks <- sapply(players, function(distr) calculate_skill(distr, players))
  ranks <- ranks[order(ranks, decreasing = TRUE)]
  print(ranks)
  
  graph_data <- data.frame(
    skill = ranks,
    player = names(ranks)
  )
  
  players_IC <- sapply(players, function(distr) {
    res <- inverse_cdf(distr)[c(11, 501, 991)] #1% and 99%
    names(res) <- c("low", "med", "up")
    res
  })
  
  graph_data2 <- data.frame(
    t(players_IC),
    player = colnames(players_IC),
    credibl = sapply(players, compute_credibility)
  )
  
  ggplot()+
    #geom_vline(aes(xintercept = skill, col = player), data=graph_data, linewidth=2)+
    geom_point(aes(x=med/100, y=credibl, col = player), data = graph_data2)+
    geom_errorbar(aes(xmin = low/100, xmax = up/100, y = credibl, col = player), data=graph_data2,
                  linewidth = 1, alpha = 0.8)+
    geom_text(aes(x=-5, y=credibl, label=player), data=graph_data2)+
    xlab("Skill")+ylab("Crédibilité")+
    theme_bw()+
    xlim(-0.1, 1)+ylim(-0.1, 1.1)
}

show_distr <- function(distr) {
  ggplot()+
    geom_point(aes(x=mu / 100, y=p), data=as.data.frame(distr))+
    xlab("Skill")+ylab("Likelihood")+
    theme_bw()+
    xlim(0, 1)
}


show_ranking_history <- function(scores) {
  
  name <- unique(unlist(scores[, c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2")]))
  name <- name[!is.na(name)]
  
  players <- list()
  for(n in name) players <- add_player(n, players)
  
  game_dates <- as.Date(unique(scores[, "date"]))
  game_dates <- sort(game_dates)
  
  drift_per_player <- lapply(name, function(nom) {
    as.character(seq.Date(as.Date(format(as.Date(min(scores[, "date"])), "%Y-%m-01")), as.Date(max(scores[, "date"])), by = "+1 month")[-1])
  })
  
  drift_dates <- as.character(seq.Date(as.Date(format(as.Date(min(scores[, "date"])), "%Y-%m-01")), as.Date(max(scores[, "date"])), by = "+1 month")[-1])
  all_dates <- c(game_dates, as.Date(drift_dates))
  all_dates <- unique(all_dates)
  all_dates <- as.Date(all_dates)
  all_dates <- sort(all_dates)
  
  graph_data <- data.frame(
    date = rep(all_dates, each = length(name)),
    player = rep(name, length(all_dates)),
    score = NA,
    played = FALSE
  )
  
  ranks <- NULL
  
  for(d in as.character(all_dates)) {
    print(d)
    
    player_in_ranking <- unique(unlist(scores[scores[, "date"] <= d, c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2")]))
    player_in_ranking <- player_in_ranking[!is.na(player_in_ranking)]
    players_today <- unique(unlist(scores[scores[, "date"] == d, c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2")]))
    players_today <- players_today[!is.na(players_today)]
    
    if(d %in% as.character(drift_dates)) {
      players[name[name %in% player_in_ranking]] <- lapply(players[name[name %in% player_in_ranking]], drift)
    }
    if(d %in% as.character(game_dates)) {
      players[player_in_ranking] <- update_scores(players=players[player_in_ranking], scores=scores[scores[, "date"] == d, ])
    }
    
    ranks <- show_current_ranking(players = players[player_in_ranking], scores = scores, init_theta = ranks)
    for(n in player_in_ranking) graph_data[graph_data[, "date"] == d & graph_data[, "player"] == n, "score"] <- ranks[n]
    for(n in players_today) graph_data[graph_data[, "date"] == d & graph_data[, "player"] == n, "played"] <- TRUE
    
    sapply(players, check_distr)
    
    print(show_detailed_skill_per_player(players[player_in_ranking]))
  }
  
  list(players, graph_data)
}

show_played_against_grid <- function(players, scores) {
  res <- matrix(0,
                nrow = length(players), ncol = length(players),
                dimnames = list(names(players), names(players))
  )
  for (i in 1:nrow(scores)) {
    res[scores[i, 3], scores[i, 4]] <- res[scores[i, 3], scores[i, 4]] + sum(scores[i, c("score_A", "score_B")])
    res[scores[i, 4], scores[i, 3]] <- res[scores[i, 4], scores[i, 3]] + sum(scores[i, c("score_A", "score_B")])
    
    # double
    if (!is.na(scores[i, 2])) {
      res[scores[i, 2], scores[i, 4]] <- res[scores[i, 2], scores[i, 4]] + sum(scores[i, c("score_A", "score_B")])
      res[scores[i, 2], scores[i, 5]] <- res[scores[i, 2], scores[i, 5]] + sum(scores[i, c("score_A", "score_B")])
      res[scores[i, 3], scores[i, 5]] <- res[scores[i, 3], scores[i, 5]] + sum(scores[i, c("score_A", "score_B")])
      res[scores[i, 4], scores[i, 2]] <- res[scores[i, 4], scores[i, 2]] + sum(scores[i, c("score_A", "score_B")])
      res[scores[i, 5], scores[i, 2]] <- res[scores[i, 5], scores[i, 2]] + sum(scores[i, c("score_A", "score_B")])
      res[scores[i, 5], scores[i, 3]] <- res[scores[i, 5], scores[i, 3]] + sum(scores[i, c("score_A", "score_B")])
    }
  }
  
  
  # ordonner selon le plus de joueurs rencontrés
  tri <- order(apply(res, 1, function(x) sum(x > 0)), decreasing = TRUE)
  
  reshaped_res <- reshape2::melt(res[tri, tri])
  reshaped_res$value <- as.logical(reshaped_res$value)
  ggplot(reshaped_res)+
    geom_tile(aes(x=Var1, y=Var2, fill = value))
}

show_ranking_history_exact <- function(scores) {
  
  name <- unique(unlist(scores[, c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2")]))
  name <- name[!is.na(name)]
  
  grid <- expand.grid(c())
  grid_id <- grid
  
  joint_density_init <- numeric(0)
  
  joint_density <- cbind(grid, p = joint_density_init)
  joint_density <- list(
    grid_id = grid_id,
    joint_distr = joint_density,
    domains = list()
  )
  
  game_dates <- as.Date(unique(scores[, "date"]))
  game_dates <- sort(game_dates)
  
  drift_dates <- as.character(seq.Date(as.Date(format(as.Date(min(scores[, "date"])), "%Y-%m-01")), as.Date(max(scores[, "date"])), by = "+1 month")[-1])
  all_dates <- c(game_dates, as.Date(drift_dates))
  all_dates <- unique(all_dates)
  all_dates <- as.Date(all_dates)
  all_dates <- sort(all_dates)
  
  graph_data <- data.frame(
    date = rep(all_dates, each = length(name)),
    player = rep(name, length(all_dates)),
    score = NA,
    played = FALSE
  )
  
  ranks <- NULL
  marginales <- list()
  nrow_before <- 0
  
  for(d in as.character(all_dates)) {
    print(d)
    
    if(nrow(joint_density$joint_distr) > 10000) {
      tmp <- simplifier_joint(joint_density, joint_density_init)
      joint_density <- tmp[[1]]
      joint_density_init <- tmp[[2]]
    }
    
    if(nrow_before != nrow(joint_density$joint_distr)) {
      print(paste0("nrow de la densité conjointe : ", nrow(
        joint_density$joint_distr), collapse = ""))
      nrow_before <- nrow(joint_density$joint_distr)
    }
    
    if(d %in% as.character(drift_dates)) {
      joint_density <- drift_exact(joint_density, joint_density_init)
    }
    if(d %in% as.character(game_dates)) {
      n_to_update <- nrow(scores[scores[, "date"] == d, ])
      i <- 0
      
      while(nrow(joint_density$joint_distr) > 10000 & i < n_to_update) {
        i <- i + 1
        
        #ajouter les nouveaux joueurs de cette partie i
        players_this_game <- unique(unlist(scores[scores[, "date"] == d, ][i, c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2"), drop = FALSE]))
        players_this_game <- players_this_game[!is.na(players_this_game)]
        
        for(n in players_this_game[!players_this_game %in% names(joint_density$domains)]) {
          print(paste0("Ajout de : ", n, collapse = ""))
          tmp <- add_player_exact(n, joint_density, joint_density_init)
          joint_density <- tmp[[1]]
          joint_density_init <- tmp[[2]]
          if(nrow_before != nrow(joint_density$joint_distr)) {
            print(paste0("nrow de la densité conjointe : ", nrow(
              joint_density$joint_distr), collapse = ""))
            nrow_before <- nrow(joint_density$joint_distr)
          }
        }
        
        joint_density <- update_scores_exact(joint_density, scores=scores[scores[, "date"] == d, ][i, , drop = FALSE])
        
        tmp <- simplifier_joint(joint_density, joint_density_init)
        joint_density <- tmp[[1]]
        joint_density_init <- tmp[[2]]
        
        if(nrow_before != nrow(joint_density$joint_distr)) {
          print(paste0("nrow de la densité conjointe : ", nrow(
            joint_density$joint_distr), collapse = ""))
          nrow_before <- nrow(joint_density$joint_distr)
        }
      }
      
      player_in_ranking <- unique(unlist(scores[scores[, "date"] <= d, c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2")]))
      player_in_ranking <- player_in_ranking[!is.na(player_in_ranking)]
      
      
      #ajouter les nouveaux joueurs
      for(n in player_in_ranking[!player_in_ranking %in% names(joint_density$domains)]) {
        tmp <- add_player_exact(n, joint_density, joint_density_init)
        joint_density <- tmp[[1]]
        joint_density_init <- tmp[[2]]
        if(nrow_before != nrow(joint_density$joint_distr)) {
          print(paste0("nrow de la densité conjointe : ", nrow(
            joint_density$joint_distr), collapse = ""))
          nrow_before <- nrow(joint_density$joint_distr)
        }
      }
      
      if(i < n_to_update) joint_density <- update_scores_exact(joint_density, scores=scores[scores[, "date"] == d, ][(i + 1):n_to_update, , drop = FALSE])
    }
    
    marginales <- marginal_from_joint(joint_density)
    
    ranks <- show_current_ranking(players = marginales[player_in_ranking], scores = scores, init_theta = ranks)
    for(n in player_in_ranking) graph_data[graph_data[, "date"] == d & graph_data[, "player"] == n, "score"] <- ranks[n]
    
    players_today <- unique(unlist(scores[scores[, "date"] == d, c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2")]))
    players_today <- players_today[!is.na(players_today)]
    
    for(n in players_today) graph_data[graph_data[, "date"] == d & graph_data[, "player"] == n, "played"] <- TRUE
    
    #print(show_detailed_skill_per_player(marginales[player_in_ranking]))
  }
  
  list(marginales, graph_data)
}

show_ranking_history_dependancy <- function(scores, dataset = "ping") {
  
  name <- unique(unlist(scores[, c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2")]))
  name <- name[!is.na(name)]
  
  clusters <- list()
  
  game_dates <- as.Date(unique(scores[, "date"]))
  game_dates <- sort(game_dates)
  
  drift_dates <- as.character(seq.Date(as.Date(format(as.Date(min(scores[, "date"])), "%Y-%m-01")), as.Date(max(scores[, "date"])), by = "+1 month")[-1])
  all_dates <- c(game_dates, as.Date(drift_dates))
  all_dates <- unique(all_dates)
  all_dates <- as.Date(all_dates)
  all_dates <- sort(all_dates)
  
  graph_data <- data.frame(
    date = rep(all_dates, each = length(name)),
    day_i = rep(order(all_dates), each = length(name)),
    player = rep(name, length(all_dates)),
    score = NA,
    rank = NA,
    played = FALSE
  )
  
  scores_players <- NULL
  marginales <- list()
  for(d in as.character(all_dates)) {
    print(d)
    
    if(d %in% as.character(drift_dates)) {
      clusters <- mapply(drift_exact2, clusters, clusters = list(clusters), SIMPLIFY = FALSE)
    }
    if(d %in% as.character(game_dates)) {
      n_to_update <- nrow(scores[scores[, "date"] == d, ])
      
      # commencer avec le simple
      # énumérer les paires de simples jouées dans la journée
      scores_subset <- scores[scores[, "date"] == d, ]
      scores_subset <- scores_subset[is.na(scores_subset$joueur_A1), ]
      if (nrow(scores_subset) > 0) {
        pairs <- mapply(function(a, b) {
          sort(c(a, b))
          }, scores_subset$joueur_A2, scores_subset$joueur_B1, SIMPLIFY = FALSE)
        unique_pairs <- unique(pairs)
        unique_pairs_game_i <- lapply(
          unique_pairs,
          function(x) which(sapply(pairs, function(y) isTRUE(all.equal(y, x))))
        )
        
        # pour chq paire, update
        for (p in seq_along(unique_pairs)) {
          # ajouter les nouveaux joueurs de cette partie i
          players_this_game <- unique(unlist(scores_subset[unique_pairs_game_i[[p]], c("joueur_A2", "joueur_B1"), drop = FALSE]))
          
          for(n in players_this_game[!players_this_game %in% unlist(sapply(clusters, `[[`, "names"))]) {
            print(paste0("Ajout de : ", n, collapse = ""))
            clusters <- add_player_dependancy(n, clusters)
          }
          
          # créer la distribution conjointe nécessaire pour la paire
          groups_to_join <- which(sapply(lapply(clusters, `[[`, "names"), function(noms) any(unique_pairs[[p]] %in% noms)))
          joint_distr_from_clusters <- join_clusters(clusters, groups_to_join)
          
          if (isFALSE(all.equal(sum(joint_distr_from_clusters$joint_distr$p), 1))) {
            print(d)
            stop("Erreur de probs")
          }
          
          # # simplifier
          # joint_distr_from_clusters <- simplifier_joint_dependancy(
          #   joint_distr_from_clusters, seuil = 0.001,
          #   absolute_max_dim = 1000000,
          #   min_no_simplif = 5000,
          #   verbose = TRUE
          # )
          
          # update la distribution
          joint_density <- update_scores_exact(
            joint_distr_from_clusters,
            scores=scores_subset[unique_pairs_game_i[[p]], , drop = FALSE]
          )
          
          # re-simplifier
          joint_density <- simplifier_joint_dependancy(
            joint_density, seuil = 0.001,
            absolute_max_dim = 200000,
            min_no_simplif = 1000,
            verbose = TRUE
          )
          
          # reclusterer
          new_clusters <- recluster_dependancy(joint_density)
          clusters <- c(new_clusters, clusters[!seq_along(clusters) %in% groups_to_join])
          print(paste0(length(clusters), " clusters, length ", paste0(sapply(clusters, function(x) length(x$names)), collapse = ", ")))
          #print(sapply(clusters, `[`, "names"))
        }
      }
      
      
      # continuer avec le double
      scores_subset <- scores[scores[, "date"] == d, ]
      scores_subset <- scores_subset[!is.na(scores_subset$joueur_A1), ]
      if (nrow(scores_subset) > 0) {
        pairs <- mapply(function(a1, a2, b1, b2) {
          sort(c(a1, a2, b1, b2))
        }, scores_subset$joueur_A1, scores_subset$joueur_A2,
        scores_subset$joueur_B1, scores_subset$joueur_B2, SIMPLIFY = FALSE)
        unique_pairs <- unique(pairs)
        unique_pairs_game_i <- lapply(
          unique_pairs,
          function(x) which(sapply(pairs, function(y) isTRUE(all.equal(y, x))))
        )
        
        # pour chq paire, update
        for (p in seq_along(unique_pairs)) {
          # ajouter les nouveaux joueurs de cette partie i
          players_this_game <- unique(unlist(scores_subset[unique_pairs_game_i[[p]], c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2"), drop = FALSE]))
          
          for(n in players_this_game[!players_this_game %in% unlist(sapply(clusters, `[[`, "names"))]) {
            print(paste0("Ajout de : ", n, collapse = ""))
            clusters <- add_player_dependancy(n, clusters)
          }
          
          # créer la distribution conjointe nécessaire pour la paire
          groups_to_join <- which(sapply(lapply(clusters, `[[`, "names"), function(noms) any(unique_pairs[[p]] %in% noms)))
          joint_distr_from_clusters <- join_clusters(clusters, groups_to_join)
          
          # # simplifier
          # joint_distr_from_clusters <- simplifier_joint_dependancy(
          #   joint_distr_from_clusters, seuil = 0.001,
          #   absolute_max_dim = 1000000,
          #   min_no_simplif = 5000,
          #   verbose = TRUE
          # )
          
          # update la distribution
          joint_density <- update_scores_exact(
            joint_distr_from_clusters,
            scores=scores_subset[unique_pairs_game_i[[p]], , drop = FALSE]
          )
          
          # re-simplifier
          joint_density <- simplifier_joint_dependancy(
            joint_density, seuil = 0.001,
            absolute_max_dim = 200000,
            min_no_simplif = 1000,
            verbose = TRUE
          )
          
          # reclusterer
          new_clusters <- recluster_dependancy(joint_density)
          clusters <- c(new_clusters, clusters[!seq_along(clusters) %in% groups_to_join])
          print(paste0(length(clusters), " clusters, length ", paste0(sapply(clusters, function(x) length(x$names)), collapse = ", ")))
          #print(sapply(clusters, `[`, "names"))
        }
      }
    }
    
    marginales <- marginal_from_joint_dependancy(clusters)
    
    scores_players <- show_current_ranking(players = marginales, scores = scores, init_theta = scores_players)
    for(n in names(marginales)) {
      graph_data[graph_data[, "date"] == d & graph_data[, "player"] == n, "score"] <- scores_players[n]
      graph_data[graph_data[, "date"] == d & graph_data[, "player"] == n, "rank"] <- which(names(sort(scores_players, decreasing = TRUE)) == n)
    }
    players_today <- unique(unlist(scores[scores[, "date"] == d, c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2")]))
    players_today <- players_today[!is.na(players_today)]
    
    for(n in players_today) graph_data[graph_data[, "date"] == d & graph_data[, "player"] == n, "played"] <- TRUE
    
    #print(show_detailed_skill_per_player(marginales))
  }
  
  list(marginales, graph_data)
}

show_played_against_grid <- function(players, scores) {
  res <- matrix(0,
                nrow = length(players), ncol = length(players),
                dimnames = list(names(players), names(players))
  )
  for (i in 1:nrow(scores)) {
    res[scores[i, 3], scores[i, 4]] <- res[scores[i, 3], scores[i, 4]] + sum(scores[i, c("score_A", "score_B")])
    res[scores[i, 4], scores[i, 3]] <- res[scores[i, 4], scores[i, 3]] + sum(scores[i, c("score_A", "score_B")])
    
    # double
    if (!is.na(scores[i, 2])) {
      res[scores[i, 2], scores[i, 4]] <- res[scores[i, 2], scores[i, 4]] + sum(scores[i, c("score_A", "score_B")])
      res[scores[i, 2], scores[i, 5]] <- res[scores[i, 2], scores[i, 5]] + sum(scores[i, c("score_A", "score_B")])
      res[scores[i, 3], scores[i, 5]] <- res[scores[i, 3], scores[i, 5]] + sum(scores[i, c("score_A", "score_B")])
      res[scores[i, 4], scores[i, 2]] <- res[scores[i, 4], scores[i, 2]] + sum(scores[i, c("score_A", "score_B")])
      res[scores[i, 5], scores[i, 2]] <- res[scores[i, 5], scores[i, 2]] + sum(scores[i, c("score_A", "score_B")])
      res[scores[i, 5], scores[i, 3]] <- res[scores[i, 5], scores[i, 3]] + sum(scores[i, c("score_A", "score_B")])
    }
  }
  
  
  # ordonner selon le plus de joueurs rencontrés
  tri <- order(apply(res, 1, function(x) sum(x > 0)), decreasing = TRUE)
  
  reshaped_res <- reshape2::melt(res[tri, tri])
  reshaped_res$value <- as.logical(reshaped_res$value)
  ggplot(reshaped_res)+
    geom_tile(aes(x=Var1, y=Var2, fill = value))
}
