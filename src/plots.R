library(ggplot2)
library(igraph)
library(profvis)

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

# credibilise <- function(distr, players, seuil = 0.7) {
#   z <- min(1, compute_credibility(distr) / seuil)
#   complement <- complement_credibilite(players)
#   rbind(cbind(mu = distr[, "mu"], p = z * distr[, "p"]), cbind(mu = complement[, "mu"], p = (1 - z) * complement[, "p"]))
# }

show_current_probs_exact <- function(clusters, ranks, pairs) {
  groups_to_join <- seq_along(clusters)
  joint_distr_from_clusters <- join_clusters(clusters, groups_to_join, seuil_simplif = 0)
  
  prob <- apply(pairs, 1, function(pair) {
    skills_idx <- joint_distr_from_clusters$grid_id[names(players)[pair]]
    MS1_list <- lapply(joint_distr_from_clusters$domains[[names(players)[pair[1]]]] / 100, transition_matrix)
    MS2_list <- lapply(joint_distr_from_clusters$domains[[names(players)[pair[2]]]] / 100, transition_matrix)
    
    MS1 <- MS1_list[skills_idx[[1]]]
    MS2 <- MS2_list[skills_idx[[2]]]
    
    sum(mapply(prob_win_point_1vs1_knowing_skills,
               MS1, MS2) * joint_distr_from_clusters$joint_distr$p)
  })
  data.frame(A=names(ranks)[pairs[, 1]], B=names(ranks)[pairs[, 2]],
             prob = round(prob, 3), prob_win_11 = round(sapply(prob, function(p) p_win_game_of(p, 11)), 3))
}

show_current_probs_exact2 <- function(clusters, ranks, pairs) {
  probs <- data.frame(
    A=names(ranks)[pairs[, 1]], B=names(ranks)[pairs[, 2]],
    prob = NA, prob_win_11 = NA
  )
  
  if (length(clusters) > 1) {
    cluster_pairs <- matrix(
      c(
        rep(1:length(clusters), each = length(clusters)),
        rep(1:length(clusters), length(clusters))
      ),
      nrow = 2,
      byrow = TRUE
    )
    cluster_pairs <- cluster_pairs[, cluster_pairs[1, ] <= cluster_pairs[2, ]]
  } else {
    cluster_pairs <- matrix(c(1, 1), ncol = 1)
  }
  
  tmp <- apply(cluster_pairs, 2, function(clust_pair) {
    noms_a_pairer <- lapply(clusters[unique(clust_pair)], `[[`, "names")
    if (length(unique(clust_pair)) == 1) {
      if (length(noms_a_pairer[[1]]) == 1) return(NULL) #si on a un joueur avec lui-même
      res <- as.data.frame(t(combn(noms_a_pairer[[1]], 2)))
      colnames(res) <- c("Var1", "Var2")
    } else {
      res <- expand.grid(noms_a_pairer)
    }
    
    if (prod(sapply(clusters[unique(clust_pair)], function(x) nrow(x$grid_id))) <= 100000) {
      joint_distr_from_clusters <- join_clusters(
        clusters, unique(clust_pair)
      )
      
      prob_res <- apply(res, 1, function(noms_pair) {
        skills_idx <- joint_distr_from_clusters$grid_id[noms_pair]
        MS1_list <- lapply(joint_distr_from_clusters$domains[[noms_pair[1]]] / 100, transition_matrix)
        MS2_list <- lapply(joint_distr_from_clusters$domains[[noms_pair[2]]] / 100, transition_matrix)
        
        MS1 <- MS1_list[skills_idx[[1]]]
        MS2 <- MS2_list[skills_idx[[2]]]
        
        likelihood_res <- numeric(length(joint_distr_from_clusters$joint_distr$p))
        same_force <- joint_distr_from_clusters$domains[[noms_pair[1]]][skills_idx[[1]]] == joint_distr_from_clusters$domains[[noms_pair[2]]][skills_idx[[2]]]
        likelihood_res[same_force] <- 0.5
        likelihood_res[!same_force] <- mapply(
          prob_win_point_1vs1_knowing_skills,
          MS1[!same_force], MS2[!same_force]
        )
        sum(likelihood_res * joint_distr_from_clusters$joint_distr$p)
      })
    } else {
      prob_res <- apply(res, 1, function(noms_pair) {
        
        joint_distr_from_clusters <- join_clusters_pair_marginal(clusters, names)
        
        skills_idx <- joint_distr_from_clusters$grid_id[noms_pair]
        MS1_list <- lapply(joint_distr_from_clusters$domains[[noms_pair[1]]] / 100, transition_matrix)
        MS2_list <- lapply(joint_distr_from_clusters$domains[[noms_pair[2]]] / 100, transition_matrix)
        
        likelihood_res <- numeric(length(joint_distr_from_clusters$joint_distr$p))
        same_force <- joint_distr_from_clusters$domains[[noms_pair[1]]][skills_idx[[1]]] == joint_distr_from_clusters$domains[[noms_pair[2]]][skills_idx[[2]]]
        likelihood_res[same_force] <- 0.5
        likelihood_res[!same_force] <- mapply(
          prob_win_point_1vs1_knowing_skills,
          MS1[!same_force], MS2[!same_force]
        )
        sum(likelihood_res * joint_distr_from_clusters$joint_distr$p)
      })
    }
    
    cbind(res, prob_res)
  }, simplify = FALSE)
  
  probs_df_not_sorted <- do.call(rbind, tmp)
  
  prob <- apply(probs, 1, function(x) {
    keep <- probs_df_not_sorted[, 1] == x[1] & probs_df_not_sorted[, 2] == x[2]
    if (sum(keep) == 1) {
      probs_df_not_sorted[keep, 3]
    } else if (sum(keep) > 0) stop("Erreur, il y a des doublons de paires")
    else {
      keep <- probs_df_not_sorted[, 1] == x[2] & probs_df_not_sorted[, 2] == x[1]
      if (sum(keep) == 0) stop("erreur")
      1 - probs_df_not_sorted[keep, 3]
    }
  })
  
  probs$prob <- round(prob, 3)
  probs$prob_win_11 <- round(sapply(prob, function(p) p_win_game_of(p, 11)), 3)
  
  probs
}

show_current_ranking <- function(clusters, scores, init_theta = NULL, show_credibility = FALSE) {
  players <- marginal_from_joint_dependancy(clusters)
  ranks <- sapply(players, function(distr) calculate_skill(distr, players))
  ordre <- order(ranks, decreasing = TRUE)
  players <- players[ordre]
  ranks <- ranks[ordre]
  pairs <- t(combn(1:length(ranks), 2))
  
  init_theta <- init_theta[names(players)]
  
  probs <- show_current_probs_exact2(clusters, ranks, pairs)
  
  # Trouver les composantes fortement connexes d'un graphe
  pairings <- players_pairs(scores)
  players_in <- as.vector(unlist(sapply(clusters, `[[`, "names")))
  
  if (!is.null(pairings)) {
    pairings <- pairings[
      sapply(
        pairings,
        function(x) all(
          x %in% players_in
        )
      )
    ]
    
    edges <- matrix(unlist(pairings), nrow=2)
    g <- graph(edges, directed = FALSE)
    clust <- components(g)$membership
    clust <- clust[names(players)]
    
    contraintes <- matrix(0, nrow = max(clust) * 2, ncol = length(players))
    
    for(i in 1:max(clust)) {
      contraintes[(i - 1) * 2 + 1, clust == i] <- 1
      contraintes[(i - 1) * 2 + 2, clust == i] <- -1
    }
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
  weights <- rep(1, length(weights))
  to_optim <- function(Forces) {
    estim <- 1 / (1 + 10^(-(Forces[pairs[, 1]] - Forces[pairs[, 2]]) / 20))
    sum(weights * (estim - probs[, "prob_win_11"])^2) / sum(weights)
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
  names(res$par) <- names(players)
  
  if(sqrt(res$value) > 0.02) warning(paste0("Convergence non parfaite des scores. RMSE of probs : ", 100 * round(sqrt(res$value), 3), "%", collapse = ""))
  
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
        ylim(0, 1)+
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

show_detailed_skill_per_player <- function(players, ordered = TRUE) {
  
  ranks <- sapply(players, function(distr) calculate_skill(distr, players))
  if (ordered)
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
    if (ordered) {
      mapply(
        function(x, n) data.frame(x, "player"=n), players2, names(players2),
        SIMPLIFY = FALSE
      )[names(sort(ranks, decreasing = TRUE))]
    } else {
      mapply(
        function(x, n) data.frame(x, "player"=n), players2, names(players2),
        SIMPLIFY = FALSE
      )
    }
    
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
    geom_hline(aes(yintercept = 0), data = graph_data2[!graph_data2$player %in% graph_data2_one_point$player, ])+
    xlab("Force")+ylab("Likelihood")+
    theme_bw()+
    xlim(0, 1)+
    scale_x_continuous(breaks = round(init_distr()[, "mu"] / 100, 2))+
    facet_grid(rows = vars(player), scales = "free_y")+
    expand_limits(y = 0)+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "none",
          panel.grid.minor.x = element_blank())
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

simplify_all <- function(clusters, dataset, ...) {
  lapply(clusters, function(joint_density) {
    simplifier_core(
      joint_density,
      dataset
    )
  })
}

simplifier_core <- function(joint_density, dataset, ...) {
  simplifier_joint_dependancy(
    joint_density,
    seuil = if(dataset == "ping") {
      1 - max(0.99, min(0.999, (0.99 - 0.999)/(500000-10000) * (nrow(joint_density$joint_distr) - 10000) + 0.999))
    } else if (dataset == "spike") {
      1 - max(0.99, min(0.999, (0.99 - 0.999)/(500000-10000) * (nrow(joint_density$joint_distr) - 10000) + 0.999))
    } else if(dataset == "pickle") {
      1 - max(0.99, min(0.999, (0.99 - 0.9995)/(500000-10000) * (nrow(joint_density$joint_distr) - 10000) + 0.9995))
    },
    absolute_max_dim = 1000000,
    min_no_simplif = if (dataset == "ping") {
      5000
    } else if (dataset == "spike") {
      10000
    } else if (dataset == "pickle") {
      10000
    },
    verbose = FALSE
  )
}

show_ranking_history_dependancy <- function(scores, dataset = "ping") {
  assign("probs_kept_counter", 1, envir = .GlobalEnv)
  
  name <- unique(unlist(scores[, c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2")]))
  name <- name[!is.na(name)]
  
  clusters <- list()
  excluded_members <- NULL
  
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
    skill = NA,
    rank_skill = NA,
    played = FALSE
  )
  
  scores_players <- NULL
  marginales <- list()
  # profvis({
  for(d in as.character(all_dates)) {
    print(d)
    print(paste0("probs kept : ", round(probs_kept_counter, 6), collapse = ""))
    
    lapply(clusters, function(joint_density) {
      if (any(is.na(joint_density$grid_id))) stop("NA dans grid_id")
      if (any(is.na(joint_density$joint_distr))) stop("NA dans joint_distr")
      if (any(mapply(
        function(grid_id, dom_len) max(grid_id) > dom_len,
        lapply(1:ncol(joint_density$grid_id), function(j) joint_density$grid_id[, j]),
        lapply(joint_density$domains, function(dom) length(dom))
      ))) stop("incohérence grid_id")
      if (any(table(unlist(lapply(clusters, `[[`, "names"))) > 1)) stop("joueurs en doublons dans clusters")
    })
    
    if (d %in% as.character(drift_dates)) {
      clusters <- mapply(drift_exact3, clusters, list(clusters), SIMPLIFY = FALSE)
    }
    if (d %in% as.character(game_dates)) {
      # commencer avec le simple
      # énumérer les paires de simples jouées dans la journée
      scores_subset <- scores[scores[, "date"] == d, ]
      scores_subset <- scores_subset[is.na(scores_subset$joueur_A1), ]
      if (nrow(scores_subset) > 0) {
        marginales <- marginal_from_joint_dependancy(clusters)
        if (length(marginales) == 0) {
          credibl <- numeric()
        } else credibl <- sapply(marginales, compute_credibility)
        
        pairs <- mapply(function(a, b) {
          sort(c(a, b))
        }, scores_subset$joueur_A2, scores_subset$joueur_B1, SIMPLIFY = FALSE)
        
        new_players <- unique(unlist(pairs)[!unlist(pairs) %in% names(marginales)])
        new_players_vec <- rep(0, length(new_players))
        names(new_players_vec) <- new_players
        credibl <- c(credibl, new_players_vec)
        
        unique_pairs <- unique(pairs)
        pairs_pseudo_credibl <- sapply(unique_pairs, function(noms) sum(credibl[noms]))
        unique_pairs_game_i <- lapply(
          unique_pairs,
          function(x) which(sapply(pairs, function(y) isTRUE(all.equal(y, x))))
        )
        
        # trier en ordre d'info
        # pour chq paire, update
        for (p in order(pairs_pseudo_credibl, decreasing = TRUE)) {
          # ajouter les nouveaux joueurs de cette partie i
          players_this_game <- unique(unlist(scores_subset[unique_pairs_game_i[[p]], c("joueur_A2", "joueur_B1"), drop = FALSE]))
          
          new_players_to_add <- players_this_game[!players_this_game %in% unlist(sapply(clusters, `[[`, "names"))]
          if (length(new_players_to_add) > 0) {
            do_simplify <- TRUE
            for(n in new_players_to_add) {
              if (length(clusters) > 0) {
                clusters <- simplify_all(clusters, dataset, probs_kept_counter = probs_kept_counter)
              }
              print(paste0("Ajout de : ", n, collapse = ""))
              clusters <- add_player_dependancy(n, clusters)
            }
          }
          
          # créer la distribution conjointe nécessaire pour la paire
          groups_to_join <- which(sapply(lapply(clusters, `[[`, "names"), function(noms) any(unique_pairs[[p]] %in% noms)))
          if (length(groups_to_join) > 1) {
            do_simplify <- TRUE
          }
          joint_distr_from_clusters <- join_clusters(clusters, groups_to_join)
          joint_distr_from_clusters$names <- names(joint_distr_from_clusters$domains)
          
          if (isFALSE(all.equal(sum(joint_distr_from_clusters$joint_distr$p), 1))) {
            print(d)
            stop("Erreur de probs")
          }
          
          # update la distribution
          joint_density <- update_scores_exact(
            joint_distr_from_clusters,
            scores=scores_subset[unique_pairs_game_i[[p]], , drop = FALSE],
            dataset
          )
          
          # re-simplifier
          if (do_simplify) {
            joint_density <- simplifier_core(joint_density, dataset, probs_kept_counter = probs_kept_counter)
            do_simplify <- FALSE
          } else {
            joint_density <- simplifier_joint_dependancy(
              joint_density,
              seuil = 0.0001,
              absolute_max_dim = Inf,
              verbose = FALSE
            )
          }
          
          # reclusterer
          tmp <- recluster_dependancy(joint_density, dataset, excluded_members, joint_distr_size_skip = Inf)
          new_clusters <- tmp[[1]]
          excluded_members <- tmp[[2]]
          clusters <- c(new_clusters, clusters[!seq_along(clusters) %in% groups_to_join])
          #print(sapply(clusters, `[`, "names"))
        }
        print(paste0(length(clusters), " clusters, length ", paste0(sapply(clusters, function(x) length(x$names)), collapse = ", ")))
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
          
          new_players_to_add <- players_this_game[!players_this_game %in% unlist(sapply(clusters, `[[`, "names"))]
          if (length(new_players_to_add) > 0) {
            do_simplify <- TRUE
            for(n in new_players_to_add) {
              if (length(clusters) > 0) {
                clusters <- simplify_all(clusters, dataset, probs_kept_counter = probs_kept_counter)
              }
              print(paste0("Ajout de : ", n, collapse = ""))
              clusters <- add_player_dependancy(n, clusters)
            }
          }
          
          # créer la distribution conjointe nécessaire pour la paire
          groups_to_join <- which(sapply(lapply(clusters, `[[`, "names"), function(noms) any(unique_pairs[[p]] %in% noms)))
          joint_distr_from_clusters <- join_clusters(clusters, groups_to_join)
          
          # update la distribution
          joint_density <- update_scores_exact(
            joint_distr_from_clusters,
            scores=scores_subset[unique_pairs_game_i[[p]], , drop = FALSE],
            dataset
          )
          
          # re-simplifier
          if (do_simplify) {
            joint_density <- simplifier_core(joint_density, dataset, probs_kept_counter = probs_kept_counter)
            do_simplify <- FALSE
          } else {
            joint_density <- simplifier_joint_dependancy(
              joint_density,
              seuil = 0.0001,
              absolute_max_dim = Inf,
              verbose = FALSE
            )
          }
          
          # reclusterer
          tmp <- recluster_dependancy(joint_density, dataset, excluded_members, joint_distr_size_skip = Inf)
          new_clusters <- tmp[[1]]
          excluded_members <- tmp[[2]]
          clusters <- c(new_clusters, clusters[!seq_along(clusters) %in% groups_to_join])
          #print(sapply(clusters, `[`, "names"))
        }
        print(paste0(length(clusters), " clusters, length ", paste0(sapply(clusters, function(x) length(x$names)), collapse = ", ")))
      }
    }
    
    print("update done")
    
    if (dataset == "ping") {
      clusters <- lapply(clusters, function(joint_density) {
        simplifier_joint_dependancy(
          joint_density, seuil = 0.0001,
          absolute_max_dim = 1000000,
          min_no_simplif = 1000,
          verbose = FALSE
        )
      })
    }
    
    marginales <- marginal_from_joint_dependancy(clusters)
    
    scores_players <- show_current_ranking(clusters, scores = scores, init_theta = scores_players)
    for(n in names(marginales)) {
      graph_data[graph_data[, "date"] == d & graph_data[, "player"] == n, "score"] <- scores_players[n]
      graph_data[graph_data[, "date"] == d & graph_data[, "player"] == n, "rank"] <- which(names(sort(scores_players, decreasing = TRUE)) == n)
      graph_data[graph_data[, "date"] == d & graph_data[, "player"] == n, "skill"] <- sum(marginales[[n]][, "mu"] * marginales[[n]][, "p"])
      graph_data[graph_data[, "date"] == d & graph_data[, "player"] == n, "rank_skill"] <- which(names(sort(sapply(marginales, function(distr) sum(distr[, "mu"] * distr[, "p"])), decreasing = TRUE)) == n)
    }
    players_today <- unique(unlist(scores[scores[, "date"] == d, c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2")]))
    players_today <- players_today[!is.na(players_today)]
    
    for(n in players_today) graph_data[graph_data[, "date"] == d & graph_data[, "player"] == n, "played"] <- TRUE
    
    print(show_detailed_skill_per_player(marginales))
  }
  # }, interval = 0.1)
  
  print(paste0("probs kept : ", round(probs_kept_counter, 6), collapse = ""))
  list(marginales, graph_data, clusters)
}

generate_GIF_images <- function(scores, dataset = "ping") {
  
  name <- unique(unlist(scores[, c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2")]))
  name <- name[!is.na(name)]
  
  clusters <- list()
  excluded_members <- NULL
  
  game_dates <- as.Date(unique(scores[, "date"]))
  game_dates <- sort(game_dates)
  
  drift_dates <- as.character(seq.Date(as.Date(format(as.Date(min(scores[, "date"])), "%Y-%m-01")), as.Date(max(scores[, "date"])), by = "+1 month")[-1])
  all_dates <- c(game_dates, as.Date(drift_dates))
  all_dates <- unique(all_dates)
  all_dates <- as.Date(all_dates)
  all_dates <- sort(all_dates)
  
  scores_players <- NULL
  marginales <- list()
  #profvis({
  # commencer avec le simple
  # énumérer les paires de simples jouées dans la journée
  scores_subset <- scores[apply(scores[, 3:4], 1, function(x) all(x %in% c("Vic", "Phil", "Will", "Éti"))), ]
  scores_subset <- scores_subset[is.na(scores_subset$joueur_A1), ]
  if (nrow(scores_subset) > 0) {
    marginales <- marginal_from_joint_dependancy(clusters)
    if (length(marginales) == 0) {
      credibl <- numeric()
    } else credibl <- sapply(marginales, compute_credibility)
    
    pairs <- mapply(function(a, b) {
      sort(c(a, b))
    }, scores_subset$joueur_A2, scores_subset$joueur_B1, SIMPLIFY = FALSE)
    
    new_players <- unique(unlist(pairs)[!unlist(pairs) %in% names(marginales)])
    new_players_vec <- rep(0, length(new_players))
    names(new_players_vec) <- new_players
    credibl <- c(credibl, new_players_vec)
    
    unique_pairs <- unique(pairs)
    pairs_pseudo_credibl <- sapply(unique_pairs, function(noms) sum(credibl[noms]))
    unique_pairs_game_i <- lapply(
      unique_pairs,
      function(x) which(sapply(pairs, function(y) isTRUE(all.equal(y, x))))
    )
    
    # tout de suite ajouter tout le monde
    for (p in order(pairs_pseudo_credibl, decreasing = TRUE)) {
      # ajouter les nouveaux joueurs de cette partie i
      players_this_game <- unique(unlist(scores_subset[unique_pairs_game_i[[p]], c("joueur_A2", "joueur_B1"), drop = FALSE]))
      
      new_players_to_add <- players_this_game[!players_this_game %in% unlist(sapply(clusters, `[[`, "names"))]
      if (length(new_players_to_add) > 0) {
        for(n in new_players_to_add) {
          print(paste0("Ajout de : ", n, collapse = ""))
          clusters <- add_player_dependancy(n, clusters)
        }
      }
    }
    
    marginales <- marginal_from_joint_dependancy(clusters)
    ggsave(
      "GIF_images_drift3/posteriori00000.png",
      plot = show_detailed_skill_per_player(marginales, FALSE),
      width = 8,
      height = 6,
      dpi = 300
    )
    
    groups_to_join <- which(sapply(lapply(clusters, `[[`, "names"), function(noms) any(unlist(unique_pairs) %in% noms)))
    joint_distr_from_clusters <- join_clusters(clusters, groups_to_join)
    joint_distr_from_clusters$names <- names(joint_distr_from_clusters$domains)
    
    if (isFALSE(all.equal(sum(joint_distr_from_clusters$joint_distr$p), 1))) {
      print(d)
      stop("Erreur de probs")
    }
    
    # update la distribution
    for (i in 1:nrow(scores_subset)) {
      groups_to_join <- which(sapply(lapply(clusters, `[[`, "names"), function(noms) any(unlist(unique_pairs) %in% noms)))
      joint_distr_from_clusters <- join_clusters(clusters, groups_to_join)
      joint_distr_from_clusters$names <- names(joint_distr_from_clusters$domains)
      
      joint_density <- update_scores_exact(
        joint_distr_from_clusters,
        scores = scores_subset[i, , drop = FALSE],
        dataset
      )
      
      # reclusterer
      tmp <- recluster_dependancy(joint_density, dataset, excluded_members, joint_distr_size_skip = Inf)
      new_clusters <- tmp[[1]]
      excluded_members <- tmp[[2]]
      clusters <- c(new_clusters, clusters[!seq_along(clusters) %in% groups_to_join])
      print(paste0(length(clusters), " clusters, length ", paste0(sapply(clusters, function(x) length(x$names)), collapse = ", ")))
      #print(sapply(clusters, `[`, "names"))
      marginales <- marginal_from_joint_dependancy(clusters)
      ggsave(
        paste0("GIF_images_drift3/posteriori", paste0(rep("0", 5-nchar(i)), collapse = ""), i, ".png", collapse = ""),
        plot = show_detailed_skill_per_player(marginales, FALSE),
        width = 8,
        height = 6,
        dpi = 300
      )
    }
    
    ggsave(
      "GIF_images2_drift3/posteriori00000.png",
      plot = show_detailed_skill_per_player(marginales, FALSE),
      width = 8,
      height = 6,
      dpi = 300
    )
    
    for (i in 1:24) {
      clusters <- mapply(drift_exact3, clusters, list(clusters), SIMPLIFY = FALSE)
      marginales <- marginal_from_joint_dependancy(clusters)
      ggsave(
        paste0("GIF_images2_drift3/posteriori", paste0(rep("0", 5-nchar(i)), collapse = ""), i, ".png", collapse = ""),
        plot = show_detailed_skill_per_player(marginales, FALSE),
        width = 8,
        height = 6,
        dpi = 300
      )
    }
  }
  NULL
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
