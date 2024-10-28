library(data.table)

calculate_skill <- function(distr_mu, players) {
  #distr_mu <- credibilise(distr_mu, players)
  round(sum(distr_mu[, "mu"] * distr_mu[, "p"]), 1) / 100
}

marginal_from_joint <- function(joint_density) {
  res <- lapply(names(joint_density$grid_id), function(nom) {
    idx_unique <- 1:length(joint_density$domains[[nom]])
    mu_unique <- joint_density$domains[[nom]]
    matrix(
      c(mu_unique,
        sapply(idx_unique, function(mu_idx) sum(joint_density$joint_distr[joint_density$grid_id[[nom]] == mu_idx, "p"]))
      ), ncol = 2, dimnames = list(NULL, c("mu", "p"))
    )
  })
  names(res) <- names(joint_density$grid_id)
  res
}

marginal_from_joint_dependancy <- function(clusters) {
  unlist(lapply(clusters, function(clust) {
    marginal_from_joint(clust)
  }), recursive = FALSE)
}

post_marginal_per_player <- function(posteriori) {
  if(ncol(posteriori) == 3) {
    posteriori <- setDT(as.data.frame(posteriori))
    
    #simplifier pour réduire les erreurs d'arrondi de by = .() dans data.table
    posteriori[, mu1 := round(mu1, 3)]
    posteriori[, mu2 := round(mu2, 3)]
    
    list(
      as.data.frame(posteriori[, .(p=sum(p)), by = .(mu1)][, .(mu=mu1, p=p)]),
      as.data.frame(posteriori[, .(p=sum(p)), by = .(mu2)][, .(mu=mu2, p=p)])
    )
  } else {
    posteriori <- setDT(as.data.frame(posteriori))
    
    #simplifier pour réduire les erreurs d'arrondi de by = .() dans data.table
    posteriori[, muA1 := round(muA1, 3)]
    posteriori[, muA2 := round(muA2, 3)]
    posteriori[, muB1 := round(muB1, 3)]
    posteriori[, muB2 := round(muB2, 3)]
    
    list(
      as.data.frame(posteriori[, .(p=sum(p)), by = .(muA1)][, .(mu=muA1, p=p)]),
      as.data.frame(posteriori[, .(p=sum(p)), by = .(muA2)][, .(mu=muA2, p=p)]),
      as.data.frame(posteriori[, .(p=sum(p)), by = .(muB1)][, .(mu=muB1, p=p)]),
      as.data.frame(posteriori[, .(p=sum(p)), by = .(muB2)][, .(mu=muB2, p=p)])
    )
  }
}

faster_expand.grid <- function (...) {
  nargs <- length(args <- list(...))
  cargs <- vector("list", nargs)
  iArgs <- seq_len(nargs)
  rep.fac <- 1L
  d <- lengths(args)
  orep <- prod(d)
  
  res <- matrix(NA, nrow = orep, ncol = nargs)
  
  for (i in iArgs) {
    x <- args[[i]]
    orep <- orep/d[i]
    res[, i] <- rep.int(rep.int(x, rep.int(rep.fac, d[i])), orep)
    rep.fac <- rep.fac * d[i]
  }
  res
}

p_win_game_of <- function(p, g = 7) {
  sapply(p, function(p_i) sum(dnbinom(0:(g-1), g, p_i)))
}

p_win_game_of_not_vec <- function(p, g = 7, ...) {
  sum(dnbinom(0:(g-1), g, p))
}

p_win_exact <- function(p, scoreA, scoreB, game_len, win) {
  stopifnot(all.equal(length(scoreA), length(scoreB), length(game_len), length(win)))
  if(length(p) > 1) stopifnot(length(scoreA) == 1)
  if(length(scoreA) > 1) stopifnot(length(p) == 1)
  stopifnot(all((scoreA > scoreB & win) | (scoreA < scoreB & !win)))
  
  
  sans_ecart <- pmax(scoreA, scoreB) == game_len
  
  ifelse(
    sans_ecart,
    dnbinom(scoreB, scoreA, p) * win + dnbinom(scoreA, scoreB, 1 - p) * (1 - win),
    p_win_exact_ecart(p, scoreA, scoreB, game_len)
  )
}

p_win_exact_not_vec <- function(p, scoreA, scoreB, game_len, win, serve_for_pt, ...) {
  stopifnot((scoreA > scoreB & win) | (scoreA < scoreB & !win))
  
  sans_ecart <- max(scoreA, scoreB) == game_len
  
  if(sans_ecart) return(dnbinom(scoreB, scoreA, p) * win + dnbinom(scoreA, scoreB, 1 - p) * (1 - win))
  p_win_exact_ecart(p, scoreA, scoreB, game_len)
}


p_win_exact_ecart <- function(p, scoreA, scoreB, game_len) {
  # ex 8-6 game_len 5
  # assurément on a passé par 4-4, puis ca prend 2 win consécutifs, une ou plusieurs fois
  breaks <- pmax(scoreA, scoreB) - game_len
  density_breaks1 <- (2 * p * (1 - p))^(breaks - 1) * p^2
  density_breaks2 <- (2 * p * (1 - p))^(breaks - 1) * (1 - p)^2
  density_ecart <- density_breaks1 * (scoreA > scoreB) + density_breaks2 * (scoreA < scoreB)
  dbinom(game_len - 1, (game_len - 1) * 2, p) * density_ecart
}

stopifnot(abs(p_win_exact_ecart(0.6, 3, 1, 2) - (0.6 * 0.4 * 2 * 0.6^2)) < 0.000001)
stopifnot(abs(p_win_exact_ecart(0.6, 1, 3, 2) - (0.6 * 0.4 * 2 * 0.4^2)) < 0.000001)
stopifnot(abs(p_win_exact_ecart(0.6, 4, 2, 2) - (0.6 * 0.4 * 2 * 0.6 * 0.4 * 2 * 0.6^2)) < 0.000001)
stopifnot(abs(p_win_exact_ecart(0.6, 3, 5, 2) - (0.6 * 0.4 * 2 * 0.6 * 0.4 * 2 * 0.6 * 0.4 * 2 * 0.4^2)) < 0.000001)
stopifnot(abs(p_win_exact_ecart(0.6, 5, 3, 3) - (0.6 * 0.6 * 0.4 * 0.4 * 6 * 0.6 * 0.4 * 2 * 0.6^2)) < 0.000001)

stopifnot(abs(p_win_exact(0.6, 1, 0, 1, 1) - 0.6) < 0.00001)
stopifnot(abs(p_win_exact(0.6, 0, 1, 1, 0) - 0.4) < 0.00001)
stopifnot(abs(p_win_exact(0.6, 2, 1, 2, 1) - 0.6*0.4*2*0.6) < 0.00001)
stopifnot(abs(p_win_exact(0.6, 1, 2, 2, 0) - 0.6*0.4*2*0.4) < 0.00001)
stopifnot(abs(p_win_exact(0.6, 1, 3, 3, 0) - 0.6*0.4^3*3) < 0.00001)
stopifnot(abs(p_win_exact(0.6, 2, 3, 3, 0) - 0.6^2*0.4^3 * 6) < 0.00001)
stopifnot(abs(p_win_exact(0.6, 2, 4, 3, 0) - 0.6^2*0.4^2*6*0.4^2) < 0.00001)
stopifnot(abs(p_win_exact(0.6, 4, 2, 3, 1) - 0.6^2*0.4^2*6*0.6^2) < 0.00001)
stopifnot(inherits(try(p_win_exact(0.6, 4, 2, 3, 0), silent = TRUE), "try-error")) # should be error
stopifnot(abs(p_win_exact(0.6, 5, 3, 3, 1) - 0.6^2*0.4^2*6*0.6*0.4*2*0.6^2) < 0.00001)
stopifnot(abs(p_win_exact(0.6, 3, 5, 3, 0) - 0.6^2*0.4^2*6*0.6*0.4*2*0.4^2) < 0.00001)
stopifnot(abs(p_win_exact(0.6, 1, 4, 4, 0) - 0.6*0.4^4*4) < 0.00001)


transition_matrix <- function(skill, k = 4, skill_fun = skill_to_p) {
  P <- matrix(NA, nrow = k+1, ncol = k+1, dimnames = list(paste0("reçu F", 0:k), paste0("frappé F", 0:k)))
  P[-1, ] <- t(sapply(1:k, function(ki) dbinom(0:k, k, skill_fun(skill, k, ki))))
  P[1, k+1] <- 1
  P[1, -(k+1)] <- 0
  P[k+1, 1] <- 1
  P[k+1, -1] <- 0
  P
}

skill_to_p <- function(skill, k, ki) {
  skill * ((k - ki) / (k - 1))
}

skill_to_p2 <- function(skill, k, ki) {
  ((skill <= 0.5) * (2 * skill * ((k - ki) / (k - 1))) + (skill > 0.5) * (2 * (skill - 1) * (ki - 1) / (k - 1) + 1))
}

prob_win_point_1vs1_knowing_skills <- function(MA, MB, k = 4) {
  
  n <- nrow(MA)
  
  M <- matrix(0, 2*n, 2*n)#, dimnames = list(
  # c(paste0("A ", rownames(MA)), paste0("B ", rownames(MB))),
  # c(paste0("B ", colnames(MA)), paste0("A ", colnames(MB)))
  # ))
  M[1:n, -1:-n] <- MA
  M[-1:-n, 1:n] <- MB
  
  for(i in 1:6) {
    M <- M %*% M
  }
  
  # prob que A gagne sachant que A sert + prob que A gagne sachant que B sert.
  # service = état milieu (dépend du sport)
  mid_k <- 1 + ceiling(k/2)
  # win = B frappe F0 ou A frappe F3
  0.5 * (M[mid_k, 1]+ M[(k+1)+mid_k, 2*(k+1)])
}

prob_win_point_2vs2_knowing_skills <- function(MA1, MA2, MB1, MB2, k = 4) {
  
  n <- nrow(MA1)
  
  #A1->B1->A2->B2
  M1 <- matrix(0, 4*n, 4*n)# , dimnames = list(
  # c(paste0("A1 ", rownames(MA1)), paste0("A2 ", rownames(MA2)), paste0("B1 ", rownames(MB1)), paste0("B2 ", rownames(MB2))),
  # c(paste0("B2 ", colnames(MB1)), paste0("B1 ", colnames(MB2)), paste0("A1 ", colnames(MA1)), paste0("A2 ", colnames(MA1)))
  # ))
  
  #A1->B2->A2->B1
  M2 <- matrix(0, 4*n, 4*n)# , dimnames = list(
  # dimnames(M1)[[1]],
  # c(paste0("B1 ", colnames(MB1)), paste0("B2 ", colnames(MB2)), paste0("A2 ", colnames(MA1)), paste0("A1 ", colnames(MA1)))
  # ))
  
  
  M1[1:n, n+n + 1:n] <- MA1
  M1[n + 1:n, n+n + n + 1:n] <- MA2
  M1[n + n + 1:n, n + 1:n] <- MB1
  M1[n + n + n + 1:n, 1:n] <- MB2
  
  M2[1:n, n+n+n + 1:n] <- MA1
  M2[n + 1:n, n+n + 1:n] <- MA2
  M2[n + n + 1:n, 1:n] <- MB1
  M2[n + n + n + 1:n, n + 1:n] <- MB2
  
  
  # 2^5 = 32 échanges durant le point
  for(i in 1:5) {
    M1 <- M1 %*% M1
    M2 <- M2 %*% M2
  }
  
  # prob win sachant A sert + prob win sachant B sert + prob win sachant A sert et cycle2, etc.
  # service = état milieu (dépend du sport)
  mid_k <- 1 + ceiling(k/2)
  
  # win = B frappe F0 ou A frappe F3
  
  0.125 * (M1[mid_k, 1] + M1[k+1+mid_k, k+1+1] + M1[2*(k+1)+mid_k, 3*(k+1)] + M1[3*(k+1)+mid_k, 4*(k+1)])+
    0.125 * (M2[mid_k, 1] + M2[k+1+mid_k, k+1+1] + M2[2*(k+1)+mid_k, 3*(k+1)] + M2[3*(k+1)+mid_k, 4*(k+1)])
  
}

prob_win_point_2vs2_knowing_skills_not_ping <- function(MA1, MA2, MB1, MB2, k = 4) {
  
  n <- nrow(MA1)
  
  # on ignore les passes
  # on supposes qu'on renvoie à n'importe quel joueur uniformément
  # on ne tient pas en compte qu'on garde le service
  M1 <- matrix(0, 4*n, 4*n)#, dimnames = list(
  #   c(paste0("A1 ", rownames(MA1)), paste0("A2 ", rownames(MA2)), paste0("B1 ", rownames(MB1)), paste0("B2 ", rownames(MB2))),
  #   c(paste0("B2 ", colnames(MB1)), paste0("B1 ", colnames(MB2)), paste0("A1 ", colnames(MA1)), paste0("A2 ", colnames(MA1)))
  # ))
  
  M1[1:n, n+n + 1:n] <- 0.5 * MA1
  M1[1:n, n+n+n + 1:n] <- 0.5 * MA1
  M1[n + 1:n, n+n + n + 1:n] <- 0.5 * MA2
  M1[n + 1:n, n+n + 1:n] <- 0.5 * MA2
  M1[n + n + 1:n, n + 1:n] <- 0.5 * MB1
  M1[n + n + 1:n, 1:n] <- 0.5 * MB1
  M1[n + n + n + 1:n, 1:n] <- 0.5 * MB2
  M1[n + n + n + 1:n, n + 1:n] <- 0.5 * MB2
  
  #eigen_info <- eigen(t(M1))
  #round(Re(eigen_info$vectors %*% diag(eigen_info$values) %*% solve(eigen_info$vectors)), 1)
  
  # 5 pour 2^5=32 coups max
  for(i in 1:5) {
    M1 <- M1 %*% M1
  }
  
  # service = état milieu (dépend du sport)
  mid_k <- 1 + ceiling(k/2)
  
  0.25 * (M1[mid_k, 1] + M1[k+1+mid_k, k+1+1] + M1[2*(k+1)+mid_k, 3*(k+1)] + M1[3*(k+1)+mid_k, 4*(k+1)] +
            M1[mid_k, k+1+1] + M1[k+1+mid_k, 1] + M1[2*(k+1)+mid_k, 4*(k+1)] + M1[3*(k+1)+mid_k, 3*(k+1)])
}

likelihood_1vs1_exact_prob_win_1_pt <- function(joint_density, player_names) {
  MS1 <- lapply(joint_density$domains[[player_names[1]]] / 100, transition_matrix)
  MS2 <- lapply(joint_density$domains[[player_names[2]]] / 100, transition_matrix)
  
  apply(
    joint_density$grid_id[player_names],
    1,
    function(idx) {
      prob_win_point_1vs1_knowing_skills(MS1[[idx[1]]], MS2[[idx[2]]])
    }
  )
}

likelihood_1vs1_exact <- function(joint_density, P_1_wins_pt, score, dataset) {
  if (dataset == "ping" || dataset == "spike") {
    not_vec_fun1 <- p_win_exact_not_vec
    not_vec_fun2 <- p_win_game_of_not_vec
  }
  if (dataset == "pickle") {
    if (as.logical(score["serve_for_pt"])) {
      not_vec_fun1 <- p_win_exact_not_vec_pickle
      not_vec_fun2 <- p_win_game_of_not_vec_pickle
    } else {
      not_vec_fun1 <- p_win_exact_not_vec
      not_vec_fun2 <- p_win_game_of_not_vec
    }
  }
  
  Likelihood_fun1 <- function(p_win_1_pt) {
    not_vec_fun1(p_win_1_pt, as.numeric(score["score_A"]), as.numeric(score["score_B"]), as.numeric(score["game_len"]), as.numeric(score["win"]), pickle_estim)
  }
  
  Likelihood_fun2 <- function(p_win_1_pt) {
    p <- sapply(p_win_1_pt, not_vec_fun2, g = as.numeric(score["game_len"]), pickle_estim = pickle_estim)
    win <- as.numeric(score["win"])
    p * win + (1 - p) * (1 - win)
  }
  
  marginales <- marginal_from_joint(joint_density)
  
  if(include_exact_points | any(sapply(marginales[name], is_exact_score_used_for_player))) {
    Likelihood <- Likelihood_fun1(P_1_wins_pt)
  } else {
    Likelihood <- Likelihood_fun2(P_1_wins_pt)
  }
  
  Likelihood
}

joint_distr <- list(grid_id = data.frame(A = c(1, 1, 2, 2), B = c(1, 2, 1, 2)),
                    joint_distr = data.frame(A = c(55, 55, 60, 60), B = c(50, 60, 50, 60), p = rep(1/4, 4)),
                    domains = list(A = c(55, 60), B = c(50, 60)))
stopifnot(all.equal(
  likelihood_1vs1_exact(joint_distr,
                        likelihood_1vs1_exact_prob_win_1_pt(joint_distr, c("A", "B")),
                        c(date = NA, joueur_A1 = NA, joueur_A2 = "A", joueur_B1 = "B", joueur_B2 = NA, win = 1, score_A = 11, score_B = 7, game_len = 11, serve_for_pt = TRUE),
                        dataset = "ping"),
  likelihood_1vs1_exact(joint_distr,
                        likelihood_1vs1_exact_prob_win_1_pt(joint_distr, c("B", "A")),
                        c(date = NA, joueur_A1 = NA, joueur_A2 = "B", joueur_B1 = "A", joueur_B2 = NA, win = 0, score_A = 7, score_B = 11, game_len = 11, serve_for_pt = TRUE),
                        dataset = "ping"),
  tolerance = 1e-6
))
stopifnot(all.equal(
  p_win_exact_not_vec(0.4, 11, 7, 11, 1, TRUE),
  p_win_exact_not_vec(1 - 0.4, 7, 11, 11, 0, TRUE)
))
stopifnot(all.equal(
  p_win_exact_not_vec(0.4, 11, 7, 11, 1, FALSE),
  p_win_exact_not_vec(1 - 0.4, 7, 11, 11, 0, FALSE)
))
joint_density <- list(grid_id = data.frame(A = c(1, 1, 2, 2), B = c(1, 2, 1, 2)),
                      joint_distr = data.frame(A = c(55, 55, 60, 60), B = c(50, 60, 50, 60), p = rep(1/4, 4)),
                      domains = list(A = c(55, 60), B = c(50, 60)))
stopifnot(all.equal(
  likelihood_1vs1_exact(joint_density,
                        likelihood_1vs1_exact_prob_win_1_pt(joint_distr, c("A", "B")),
                        c(date = NA, joueur_A1 = NA, joueur_A2 = "A", joueur_B1 = "B", joueur_B2 = NA, win = 1, score_A = 11, score_B = 7, game_len = 11, serve_for_pt = FALSE),
                        dataset = "ping"),
  likelihood_1vs1_exact(joint_density,
                        likelihood_1vs1_exact_prob_win_1_pt(joint_distr, c("B", "A")),
                        c(date = NA, joueur_A1 = NA, joueur_A2 = "B", joueur_B1 = "A", joueur_B2 = NA, win = 0, score_A = 7, score_B = 11, game_len = 11, serve_for_pt = FALSE),
                        dataset = "ping"),
  tolerance = 1e-6
))

likelihood_2vs2_exact_prob_win_1_pt <- function(joint_density, player_names) {
  MSA1 <- lapply(joint_density$domains[[player_names[1]]] / 100, transition_matrix)
  MSA2 <- lapply(joint_density$domains[[player_names[2]]] / 100, transition_matrix)
  MSB1 <- lapply(joint_density$domains[[player_names[3]]] / 100, transition_matrix)
  MSB2 <- lapply(joint_density$domains[[player_names[4]]] / 100, transition_matrix)
  
  if (dataset == "ping") prob_point_fun <- prob_win_point_2vs2_knowing_skills
  if (dataset == "spike" || dataset == "pickle") prob_point_fun <- prob_win_point_2vs2_knowing_skills_not_ping
  
  apply(
    joint_density$grid_id[player_names],
    1,
    function(idx) {
      prob_point_fun(
        MSA1[[idx[1]]], MSA2[[idx[2]]], MSB1[[idx[3]]], MSB2[[idx[4]]]
      )
    }
  )
}

joint_density <- list(
  grid_id = expand.grid(A=1:2, B=1:2, C=1:2, D=1:2),
  joint_distr = cbind(
    expand.grid(
      A = c(50, 60),
      B = c(55, 62),
      C = c(40, 46), 
      D = c(52, 54)
    ),
    p = c(0.1, 0.05, 0.03, 0.05, 0.02, 0.03, 0.06, 0.03, 0.04, 0.06, 0.03, 0.04, 0.18, 0.08, 0.11, 0.09)
  ),
  domains = list(
    A = c(50, 60),
    B = c(55, 62),
    C = c(40, 46),
    D = c(52, 54)
  )
)
stopifnot(all.equal(
  likelihood_2vs2_exact_prob_win_1_pt(joint_density, c("A", "B", "C", "D")),
  likelihood_2vs2_exact_prob_win_1_pt(joint_density, c("B", "A", "C", "D"))
))
stopifnot(all.equal(
  likelihood_2vs2_exact_prob_win_1_pt(joint_density, c("A", "B", "C", "D")),
  likelihood_2vs2_exact_prob_win_1_pt(joint_density, c("A", "B", "D", "C"))
))
stopifnot(all.equal(
  likelihood_2vs2_exact_prob_win_1_pt(joint_density, c("A", "B", "C", "D")),
  likelihood_2vs2_exact_prob_win_1_pt(joint_density, c("B", "A", "D", "C"))
))
stopifnot(all.equal(
  likelihood_2vs2_exact_prob_win_1_pt(joint_density, c("A", "B", "C", "D")),
  1 - likelihood_2vs2_exact_prob_win_1_pt(joint_density, c("C", "D", "A", "B")),
  tolerance = 0.001
))

likelihood_2vs2_exact <- function(joint_density, P_A_wins_pt, score, dataset) {
  if (dataset == "ping" || dataset == "spike") {
    not_vec_fun1 <- p_win_exact_not_vec
    not_vec_fun2 <- p_win_game_of_not_vec
  }
  if (dataset == "pickle") {
    if (as.logical(score["serve_for_pt"])) {
      not_vec_fun1 <- p_win_exact_not_vec_pickle
      not_vec_fun2 <- p_win_game_of_not_vec_pickle
    } else {
      not_vec_fun1 <- p_win_exact_not_vec
      not_vec_fun2 <- p_win_game_of_not_vec
    }
  }
  Likelihood_fun1 <- function(P_A_wins_pt) {
    not_vec_fun1(P_A_wins_pt, as.numeric(score["score_A"]), as.numeric(score["score_B"]), as.numeric(score["game_len"]), as.numeric(score["win"]), pickle_estim)
  }
  
  Likelihood_fun2 <- function(P_A_wins_pt) {
    p <- sapply(P_A_wins_pt, not_vec_fun2, g = as.numeric(score["game_len"]), pickle_estim)
    win <- as.numeric(score["win"])
    p * win + (1 - p) * (1 - win)
  }
  
  marginales <- marginal_from_joint(joint_density)
  
  
  #weighter les games selon le nombre de jours passé avec (0.5^(2/365))^-x
  if(include_exact_points | any(sapply(marginales[name], is_exact_score_used_for_player))) {
    Likelihood <- Likelihood_fun1(P_A_wins_pt)
  } else {
    Likelihood <- Likelihood_fun2(P_A_wins_pt)
  }
  
  Likelihood
}

players_pairs <- function(scores) {
  # scores <- scores[is.na(scores[, "joueur_A1"]), ] #juste trouver les paires en 1 vs 1
  tmp <- apply(
    scores[, 3:4], 1, function(x) {
      names(x) <- NULL
      sort(x)
    }
  )
  if (length(tmp) > 0) {
    unique(
      lapply(
        1:ncol(tmp),
        function(i) tmp[, i]
      )
    )
  } else NULL
}

update_scores_exact <- function(joint_density, scores, dataset) {
  if (!(all(is.na(scores$joueur_A1)) | all(!is.na(scores$joueur_A1)))) stop("erreur, le simple et double doivent
                                                                         être exécutés séparément")
  # simple
  if (is.na(scores$joueur_A1[1])) {
    player_names <- unique(unlist(scores[c("joueur_A2", "joueur_B1")]))
    if (length(player_names) != 2) stop("erreur, on doit n'avoir que 2 joueurs en simple à updater en même temps")
    
    P_1_wins_pt <- likelihood_1vs1_exact_prob_win_1_pt(joint_density, player_names[1:2])
  } else {
    # double
    player_names <- unique(unlist(scores[c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2")]))
    if (length(player_names) != 4) stop("erreur, on doit n'avoir que 4 joueurs en double à updater en même temps")
    
    if (any(apply(
      scores[paste0("joueur_", c("A1", "A2", "B1", "B2"))], 1,
      function(x) all(
        sort(x[1:2]) == sort(player_names[1:2])) | all(sort(x[3:4]) == sort(player_names[1:2]))))) {
      # print("Cas 1-2")
      P_A_wins_pt_V1 <- likelihood_2vs2_exact_prob_win_1_pt(joint_density, player_names[1:4])
    }
    if (any(apply(
      scores[paste0("joueur_", c("A1", "A2", "B1", "B2"))], 1,
      function(x) all(
        sort(x[1:2]) == sort(player_names[c(1, 3)])) | all(sort(x[3:4]) == sort(player_names[c(1, 3)]))))) {
      # print("Cas 1-3")
      P_A_wins_pt_V2 <- likelihood_2vs2_exact_prob_win_1_pt(joint_density, player_names[c(1, 3, 2, 4)])
    }
    if (any(apply(
      scores[paste0("joueur_", c("A1", "A2", "B1", "B2"))], 1,
      function(x) all(
        sort(x[1:2]) == sort(player_names[c(1, 4)])) | all(sort(x[3:4]) == sort(player_names[c(1, 4)]))))) {
      # print("Cas 1-4")
      P_A_wins_pt_V3 <- likelihood_2vs2_exact_prob_win_1_pt(joint_density, player_names[c(1, 4, 2, 3)])
    }
  }
  
  
  posteriori <- apply(
    matrix(
      apply(scores, 1, function(score) {
        # print(score)
        if(is.na(score["joueur_A1"])) return(
          likelihood_1vs1_exact(
            joint_density,
            if (score["joueur_A2"] == player_names[1]) P_1_wins_pt else 1 - P_1_wins_pt,
            score,
            dataset
          )
        )
        
        case <- 0
        x <- unlist(score[paste0("joueur_", c("A1", "A2", "B1", "B2"))])
        if (all(sort(x[1:2]) == sort(player_names[1:2]))) case <- 1
        if (all(sort(x[3:4]) == sort(player_names[1:2]))) case <- 2
        if (all(sort(x[1:2]) == sort(player_names[c(1, 3)]))) case <- 3
        if (all(sort(x[3:4]) == sort(player_names[c(1, 3)]))) case <- 4
        if (all(sort(x[1:2]) == sort(player_names[c(1, 4)]))) case <- 5
        if (all(sort(x[3:4]) == sort(player_names[c(1, 4)]))) case <- 6
        
        
        if (case %in% 1:2) {
          return(
            likelihood_2vs2_exact(
              joint_density,
              if (case == 1) P_A_wins_pt_V1 else 1 - P_A_wins_pt_V1,
              score,
              dataset
            )
          )
        }
        if (case %in% 3:4) {
          return(
            likelihood_2vs2_exact(
              joint_density,
              if (case == 3) P_A_wins_pt_V2 else 1 - P_A_wins_pt_V2,
              score,
              dataset
            )
          )
        }
        if (case %in% 5:6) {
          return(
            likelihood_2vs2_exact(
              joint_density,
              if (case == 5) P_A_wins_pt_V3 else 1 - P_A_wins_pt_V3,
              score,
              dataset
            )
          )
        }
        
        
      }), nrow = nrow(joint_density$joint_distr), ncol = nrow(scores)
    ),
    1, prod
  ) * joint_density$joint_distr$p
  
  posteriori <- posteriori / sum(posteriori)
  
  joint_density$joint_distr$p <- posteriori
  joint_density
}

