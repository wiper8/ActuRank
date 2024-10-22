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


distr_F1_F2_1vs1 <- function(distr_F1, distr_F2) {
  
  tmp <- faster_expand.grid(1:nrow(distr_F1), 1:nrow(distr_F2))
  distr_F1_F2 <- cbind(distr_F1[tmp[, 1], , drop = FALSE], distr_F2[tmp[, 2], , drop = FALSE])
  colnames(distr_F1_F2) <- c(t(sapply(colnames(distr_F1), function(x) paste0(x, c("1", "2")))))
  distr_F1_F2
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

distr_F1_F2_2vs2 <- function(distr_FA1, distr_FA2, distr_FB1, distr_FB2) {
  
  tmp <- faster_expand.grid(1:nrow(distr_FA1), 1:nrow(distr_FA2), 1:nrow(distr_FB1), 1:nrow(distr_FB2))
  distr_F1_F2 <- cbind(distr_FA1[tmp[, 1], ], distr_FA2[tmp[, 2], ], distr_FB1[tmp[, 3], ], distr_FB2[tmp[, 4], ])
  colnames(distr_F1_F2) <- c(t(sapply(colnames(distr_FA1), function(x) paste0(x, c("A1", "A2", "B1", "B2")))))
  distr_F1_F2
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


distr_P_1vs1 <- function(distr_F1_F2) {
  cbind(
    distr_F1_F2,
    "P" = distr_F1_F2[, "F1"] / (distr_F1_F2[, "F1"] + distr_F1_F2[, "F2"])
  )
}

#TODO réfléchir a si c'est vraiment le calcul que je veux
#présentement c'est l'équivalent à dire chq point sont en 1vs1, mais on alternes les possibilitées
distr_P_2vs2 <- function(distr_F1_F2) {
  cbind(distr_F1_F2, 
        "P" = (distr_F1_F2[, "FA1"] / (distr_F1_F2[, "FA1"] + distr_F1_F2[, "FB1"]) + distr_F1_F2[, "FA2"] / (distr_F1_F2[, "FA2"] + distr_F1_F2[, "FB1"]) +
     distr_F1_F2[, "FA1"] / (distr_F1_F2[, "FA1"] + distr_F1_F2[, "FB2"]) + distr_F1_F2[, "FA2"] / (distr_F1_F2[, "FA2"] + distr_F1_F2[, "FB2"])) / 4
  )
}

transition_matrix <- function(skill, k = 3) {
  P <- matrix(NA, nrow = k+1, ncol = k+1, dimnames = list(paste0("reçu F", 0:k), paste0("frappé F", 0:k)))
  P[-1, ] <- t(sapply(1:k, function(ki) dbinom(0:k, k, skill * ((1 - ki) / (k - 1) + 1))))
  P[1, k+1] <- 1
  P[1, -(k+1)] <- 0
  P
}

prob_win_point_1vs1_knowing_skills <- function(MA, MB, k = 3) {
  
  M <- matrix(0, 2*nrow(MA), 2*ncol(MA))# , dimnames = list(
    # c(paste0("A ", rownames(MA)), paste0("B ", rownames(MB))),
    # c(paste0("B ", colnames(MA)), paste0("A ", colnames(MB)))
  # ))
  M[1:nrow(MA), -1:-nrow(MA)] <- MA
  M[-1:-nrow(MA), 1:nrow(MA)] <- MB
  
  for(i in 1:6) {
    M <- M %*% M
  }
  
  # prob que A gagne sachant que A sert +  prob que A gagne sachant que B sert.
  # service = état milieu (dépend du sport)
  # win = B frappe F0 ou A frappe F3
  # TODO enlever le 3 hardcodé et faire dépendre selon k
  # devrait etre état milieu environ
  
  0.5 * (M[3, 1]+ M[(k+1)+3, 2*(k+1)])
}

prob_win_point_2vs2_knowing_skills <- function(MA1, MA2, MB1, MB2, k = 3) {
  
  #A1->B1->A2->B2
  M1 <- matrix(0, 4*nrow(MA1), 4*ncol(MA1))# , dimnames = list(
    # c(paste0("A1 ", rownames(MA1)), paste0("A2 ", rownames(MA2)), paste0("B1 ", rownames(MB1)), paste0("B2 ", rownames(MB2))),
    # c(paste0("B2 ", colnames(MB1)), paste0("B1 ", colnames(MB2)), paste0("A1 ", colnames(MA1)), paste0("A2 ", colnames(MA1)))
  # ))
  
  #A1->B2->A2->B1
  M2 <- matrix(0, 4*nrow(MA1), 4*ncol(MA1))# , dimnames = list(
    # dimnames(M1)[[1]],
    # c(paste0("B1 ", colnames(MB1)), paste0("B2 ", colnames(MB2)), paste0("A2 ", colnames(MA1)), paste0("A1 ", colnames(MA1)))
  # ))
  
  
  M1[1:nrow(MA1), ncol(MB1)+ncol(MB2) + 1:ncol(MA1)] <- MA1
  M1[nrow(MA1) + 1:nrow(MA2), ncol(MB1)+ncol(MB2) + ncol(MA1) + 1:ncol(MA2)] <- MA2
  M1[nrow(MA1) + nrow(MA2) + 1:nrow(MB1), ncol(MB2) + 1:ncol(MB1)] <- MB1
  M1[nrow(MA1) + nrow(MA2) + nrow(MB2) + 1:nrow(MB2), 1:ncol(MB2)] <- MB2
  
  M2[1:nrow(MA1), ncol(MB1)+ncol(MB2)+ncol(MA2) + 1:ncol(MA1)] <- MA1
  M2[nrow(MA1) + 1:nrow(MA2), ncol(MB1)+ncol(MB2) + 1:ncol(MA2)] <- MA2
  M2[nrow(MA1) + nrow(MA2) + 1:nrow(MB1), 1:ncol(MB1)] <- MB1
  M2[nrow(MA1) + nrow(MA2) + nrow(MB2) + 1:nrow(MB2), ncol(MB1) + 1:ncol(MB2)] <- MB2
  
  #eigen_info <- eigen(t(M1))
  #round(Re(eigen_info$vectors %*% diag(eigen_info$values) %*% solve(eigen_info$vectors)), 1)
  
  # 2^5 = 32 échanges durant le point
  for(i in 1:5) {
    M1 <- M1 %*% M1
    M2 <- M2 %*% M2
  }
  
  # prob win sachant A sert + prob win sachant B sert + prob win sachant A sert et cycle2, etc.
  # service = départ dans l'état F2
  # win = B frappe F0 ou A frappe F3
  
  0.125 * (M1[3, 1] + M1[k+1+3, k+1+1] + M1[2*(k+1)+3, 3*(k+1)] + M1[3*(k+1)+3, 4*(k+1)])+
    0.125 * (M2[3, 1] + M1[k+1+3, k+1+1] + M1[2*(k+1)+3, 3*(k+1)] + M1[3*(k+1)+3, 4*(k+1)])
  
}

prob_win_point_2vs2_knowing_skills_spikeball <- function(MA1, MA2, MB1, MB2, k = 3) {
  # on ignore les passes
  # on supposes qu'on renvoie à n'importe quel joueur uniformément
  # on ne tient pas en compte qu'on garde le service
  M1 <- matrix(0, 4*nrow(MA1), 4*ncol(MA1), dimnames = list(
    c(paste0("A1 ", rownames(MA1)), paste0("A2 ", rownames(MA2)), paste0("B1 ", rownames(MB1)), paste0("B2 ", rownames(MB2))),
    c(paste0("B2 ", colnames(MB1)), paste0("B1 ", colnames(MB2)), paste0("A1 ", colnames(MA1)), paste0("A2 ", colnames(MA1)))
  ))
  
  M1[1:nrow(MA1), ncol(MB1)+ncol(MB2) + 1:ncol(MA1)] <- 0.5 * MA1
  M1[1:nrow(MA1), ncol(MB1)+ncol(MB2)+ncol(MA2) + 1:ncol(MA1)] <- 0.5 * MA1
  M1[nrow(MA1) + 1:nrow(MA2), ncol(MB1)+ncol(MB2) + ncol(MA1) + 1:ncol(MA2)] <- 0.5 * MA2
  M1[nrow(MA1) + 1:nrow(MA2), ncol(MB1)+ncol(MB2) + 1:ncol(MA2)] <- 0.5 * MA2
  M1[nrow(MA1) + nrow(MA2) + 1:nrow(MB1), ncol(MB2) + 1:ncol(MB1)] <- 0.5 * MB1
  M1[nrow(MA1) + nrow(MA2) + 1:nrow(MB1), 1:ncol(MB1)] <- 0.5 * MB1
  M1[nrow(MA1) + nrow(MA2) + nrow(MB2) + 1:nrow(MB2), 1:ncol(MB2)] <- 0.5 * MB2
  M1[nrow(MA1) + nrow(MA2) + nrow(MB2) + 1:nrow(MB2), ncol(MB1) + 1:ncol(MB2)] <- 0.5 * MB2
  
  #eigen_info <- eigen(t(M1))
  #round(Re(eigen_info$vectors %*% diag(eigen_info$values) %*% solve(eigen_info$vectors)), 1)
  
  # 5 pour 2^5=32 coups max
  for(i in 1:5) {
    M1 <- M1 %*% M1
  }
  
  0.25 * (M1[3, 1] + M1[k+1+3, k+1+1] + M1[2*(k+1)+3, 3*(k+1)] + M1[3*(k+1)+3, 4*(k+1)] +
           M1[3, k+1+1] + M1[k+1+3, 1] + M1[2*(k+1)+3, 4*(k+1)] + M1[3*(k+1)+3, 3*(k+1)])
}


posteriori_1vs1_vectorized <- function(distr_S1, distr_S2, game_len, win, date, scoreA, scoreB, name) {
  
  distr_S1_S2 <- distr_F1_F2_1vs1(distr_S1, distr_S2)
  
  MS1 <- rep(lapply(distr_S1[, "mu"] / 100, transition_matrix), nrow(distr_S2))
  MS2 <- rep(lapply(distr_S2[, "mu"] / 100, transition_matrix), each = nrow(distr_S1))
  
  
  distr_P <- cbind(mu1=distr_S1_S2[, "mu1"],
                   mu2=distr_S1_S2[, "mu2"],
                   "P_1_wins_pt"=mapply(function(MS1, MS2) prob_win_point_1vs1_knowing_skills(MS1, MS2),
                                        MS1, MS2),
                   "p_s1_s2" = distr_S1_S2[, "p1"] * distr_S1_S2[, "p2"]
  )
  
  Likelihood_fun1 <- function(p_win_1_pt) {
    p <- p_win_exact(p_win_1_pt, scoreA, scoreB, game_len, win)
    prod(p)#^((0.5^(2 / 365))^as.numeric(Sys.Date() - date)))
  }
  
  Likelihood_fun2 <- function(p_win_1_pt) {
    p <- sapply(game_len, p_win_game_of_not_vec, p=p_win_1_pt)
    
    prod((p * win + (1 - p) * (1 - win)))#^((0.5^(2 / 365))^as.numeric(Sys.Date() - date)))
  }
  
  #weighter les games selon le nombre de jours passé avec (0.5^(2/365))^-x
  if(include_exact_points | any(sapply(list(distr_S1, distr_S2), is_exact_score_used_for_player))) {
    Likelihood <- sapply(distr_P[, "P_1_wins_pt"], Likelihood_fun1) * distr_P[, "p_s1_s2"]
  } else {
    Likelihood <- sapply(distr_P[, "P_1_wins_pt"], Likelihood_fun2) * distr_P[, "p_s1_s2"]
  }
  
  Likelihood <- Likelihood / sum(Likelihood)
  
  posteriori <- cbind(
    distr_P[, c("mu1", "mu2")],
    "p"=Likelihood
  )
  
  posteriori
}

likelihood_1vs1_exact <- function(joint_density, score, dataset) {
  name <- unlist(score[c("joueur_A2", "joueur_B1")])
  
  MS1 <- lapply(joint_density$domains[[name[1]]] / 100, transition_matrix)
  MS2 <- lapply(joint_density$domains[[name[2]]] / 100, transition_matrix)
  
  P_1_wins_pt <- list(
    matrix(
      apply(
        expand.grid(1:length(MS1), 1:length(MS2)),
        1,
        function(idx) prob_win_point_1vs1_knowing_skills(MS1[[idx[1]]], MS2[[idx[2]]])
      ),
      nrow=length(MS1)
    )
  )
  #P_1_wins_pt <- list(t(sapply(1:length(MS1), function(i) sapply(1:length(MS2), function(j) prob_win_point_1vs1_knowing_skills(MS1[[i]], MS2[[j]])))))
  P_1_wins_pt <- apply(joint_density$grid_id[name], 1, function(idx) do.call(`[`, c(P_1_wins_pt, as.list(idx))))
  
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
  
  
  #weighter les games selon le nombre de jours passé avec (0.5^(2/365))^-x
  if(include_exact_points | any(sapply(marginales[name], is_exact_score_used_for_player))) {
    Likelihood <- Likelihood_fun1(P_1_wins_pt)
  } else {
    Likelihood <- Likelihood_fun2(P_1_wins_pt)
  }
  
  Likelihood
}

likelihood_2vs2_exact <- function(joint_density, score, dataset) {
  name <- unlist(score[c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2")])
  
  MSA1 <- lapply(joint_density$domains[[name[1]]] / 100, transition_matrix)
  MSA2 <- lapply(joint_density$domains[[name[2]]] / 100, transition_matrix)
  MSB1 <- lapply(joint_density$domains[[name[3]]] / 100, transition_matrix)
  MSB2 <- lapply(joint_density$domains[[name[4]]] / 100, transition_matrix)
  
  if (dataset == "ping") prob_point_fun <- prob_win_point_2vs2_knowing_skills
  if (dataset == "spike" || dataset == "pickle") prob_point_fun <- prob_win_point_2vs2_knowing_skills_spikeball
  
  P_A_wins_pt <- list(
    array(
      apply(
        expand.grid(1:length(MSA1), 1:length(MSA2), 1:length(MSB1), 1:length(MSB2)),
        1,
        function(idx) prob_point_fun(
          MSA1[[idx[1]]], MSA2[[idx[2]]], MSB1[[idx[3]]], MSB2[[idx[4]]])
      ),
      dim = c(length(MSA1), length(MSA2), length(MSB1), length(MSB2))
    )
  )
  
  P_A_wins_pt <- apply(joint_density$grid_id[name], 1, function(idx) do.call(`[`, c(P_A_wins_pt, as.list(idx))))
  
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

posteriori_1vs1 <- function(distr_S1, distr_S2, game_len, win, date, scoreA, scoreB, name) {
  
  distr_S1_S2 <- distr_F1_F2_1vs1(distr_S1, distr_S2)
  
  MS1 <- rep(lapply(distr_S1[, "mu"] / 100, transition_matrix), nrow(distr_S2))
  MS2 <- rep(lapply(distr_S2[, "mu"] / 100, transition_matrix), each = nrow(distr_S1))
  
  
  distr_P <- cbind(mu1=distr_S1_S2[, "mu1"],
                   mu2=distr_S1_S2[, "mu2"],
                   "P_1_wins_pt"=mapply(function(MS1, MS2) prob_win_point_1vs1_knowing_skills(MS1, MS2),
         MS1, MS2),
         "p_s1_s2" = distr_S1_S2[, "p1"] * distr_S1_S2[, "p2"]
  )
  
  if(include_exact_points | any(sapply(list(distr_S1, distr_S2), is_exact_score_used_for_player))) {
    distr_P <- cbind(distr_P, "P_win_game" = p_win_exact(distr_P[, "P_1_wins_pt"], scoreA, scoreB, game_len, win))
    
    #weighter les games selon le nombre de jours passé avec (0.5^(2/365))^-x
    Likelihood <- distr_P[, "P_win_game"]#^((0.5^(2/365))^as.numeric(Sys.Date() - date)) * distr_P[, "p_s1_s2"]
  } else {
    distr_P <- cbind(distr_P, "P_win_game" = p_win_game_of(distr_P[, "P_1_wins_pt"], game_len))
    
    #weighter les games selon le nombre de jours passé avec (0.5^(2/365))^-x
    Likelihood <- (distr_P[, "P_win_game"] * win + (1-distr_P[, "P_win_game"]) * (1-win))#^((0.5^(2/365))^as.numeric(Sys.Date() - date)) * distr_P[, "p_s1_s2"]
  }
  
  Likelihood <- Likelihood/sum(Likelihood)
  
  posteriori <- cbind(
    distr_P[, c("mu1", "mu2")],
    "p"=Likelihood
  )
  
  posteriori
}


posteriori_2vs2 <- function(distr_SA1, distr_SA2,
                            distr_SB1, distr_SB2,
                            game_len, win, date, scoreA, scoreB) {
  
  distr_SA1_SA2_SB1_SB2 <- distr_F1_F2_2vs2(distr_SA1, distr_SA2,
                                  distr_SB1, distr_SB2)
  
  MSA1 <- rep(lapply(distr_SA1[, "mu"] / 100, transition_matrix),
              nrow(distr_SA2) * nrow(distr_SB1) * nrow(distr_SB2))
  MSA2 <- rep(rep(lapply(distr_SA2[, "mu"] / 100, transition_matrix),
                  each = nrow(distr_SA1)), nrow(distr_SB1) * nrow(distr_SB2))
  MSB1 <- rep(rep(lapply(distr_SB1[, "mu"] / 100, transition_matrix),
                  each = nrow(distr_SA1) * nrow(distr_SA2)), nrow(distr_SB2))
  MSB2 <- rep(lapply(distr_SB2[, "mu"] / 100, transition_matrix),
                  each = nrow(distr_SA1) * nrow(distr_SA2) * nrow(distr_SB1))
  
  distr_P <- cbind(muA1=distr_SA1_SA2_SB1_SB2[, "muA1"],
                   muA2=distr_SA1_SA2_SB1_SB2[, "muA2"],
                   muB1=distr_SA1_SA2_SB1_SB2[, "muB1"],
                   muB2=distr_SA1_SA2_SB1_SB2[, "muB2"],
                   "P_A_wins_pt"=mapply(function(MSA1, MSA2, MSB1, MSB2) prob_win_point_2vs2_knowing_skills(MSA1, MSA2, MSB1, MSB2),
                                        MSA1, MSA2, MSB1, MSB2),
                   "p_sa1_sa2_sb1_sb2" = apply(distr_SA1_SA2_SB1_SB2[, c("pA1", "pA2", "pB1", "pB2")], 1, prod)
  )
  
  if(include_exact_points | any(sapply(list(distr_SA1, distr_SA2, distr_SB1, distr_SB2), is_exact_score_used_for_player))) {
    distr_P <- cbind(distr_P, "P_win_game" = p_win_exact(distr_P[, "P_A_wins_pt"], scoreA, scoreB, game_len, win))
    #weighter les games selon le nombre de jours passé avec (0.5^(2/365))^-x
    Likelihood <- distr_P[, "P_win_game"]#^((0.5^(2/365))^as.numeric(Sys.Date() - date)) * distr_P[, "p_sa1_sa2_sb1_sb2"]
  } else {
    distr_P <- cbind(distr_P, "P_win_game" = p_win_game_of(distr_P[, "P_A_wins_pt"], game_len))
    #weighter les games selon le nombre de jours passé avec (0.5^(2/365))^-x
    Likelihood <- (distr_P[, "P_win_game"] * win + (1-distr_P[, "P_win_game"]) * (1-win))#^((0.5^(2/365))^as.numeric(Sys.Date() - date)) * distr_P[, "p_sa1_sa2_sb1_sb2"]
  }
  
  Likelihood <- Likelihood/sum(Likelihood)
  Likelihood <- cbind(Likelihood, distr_P)
  
  posteriori <- cbind(
    Likelihood[, c("muA1", "muA2", "muB1", "muB2")],
    "p"=Likelihood[, "Likelihood"]
  )
  
  posteriori
}

posteriori_of_game_simplified_vectorized <- function(players, scores) {
  A <- names(players)[1]; B <- names(players)[2]
  
  players[[A]] <- round_up_domain(players[[A]])
  players[[B]] <- round_up_domain(players[[B]])
  
  
  posteriori <- posteriori_1vs1_vectorized(
    distr_S1 = players[[A]],
    distr_S2 = players[[B]],
    game_len = as.numeric(scores[, "game_len"]),
    win = as.numeric(scores[, "win"]) * (scores[, "joueur_A2"] == A) + (1 - as.numeric(scores[, "win"])) * (scores[, "joueur_A2"] == B),
    date = as.Date(scores[, "date"]),
    scoreA = as.numeric(scores[, "score_A"]) * (scores[, "joueur_A2"] == A) + as.numeric(scores[, "score_B"]) * (scores[, "joueur_A2"] == B),
    scoreB = as.numeric(scores[, "score_B"]) * (scores[, "joueur_A2"] == A) + as.numeric(scores[, "score_A"]) * (scores[, "joueur_A2"] == B),
    name = c(A, B)
  )
  posteriori_per_player <- post_marginal_per_player(posteriori)
  
  players[[A]] <- as.matrix(posteriori_per_player[[1]])
  players[[B]] <- as.matrix(posteriori_per_player[[2]])
  players
}

posteriori_of_game_simplified <- function(players, score) {
  if(is.na(score[, "joueur_A1"])) {
    
    players[[score[, "joueur_A2"]]] <- round_up_domain(players[[score[, "joueur_A2"]]])
    players[[score[, "joueur_B1"]]] <- round_up_domain(players[[score[, "joueur_B1"]]])
    
    tmp <- tails_simplifier_1vs1(distr1 = players[[score[, "joueur_A2"]]],
                                 distr2 = players[[score[, "joueur_B1"]]])
    
    probs_ignorees1 <- sum(players[[score[, "joueur_A2"]]][!tmp[["keep1"]], "p"])
    probs_ignorees2 <- sum(players[[score[, "joueur_B1"]]][!tmp[["keep2"]], "p"])
    
    if(probs_ignorees1 > 0.01 | probs_ignorees2 > 0.01) {
      probs_ignorees1 <- 0
      probs_ignorees2 <- 0
      tmp <- list(keep1 = rep(T, length(tmp[["keep1"]])),
                  keep2 = rep(T, length(tmp[["keep2"]])))
    } else {
      #print(paste0("  ", prod(sapply(tmp, sum))))
    }
    
    posteriori <- posteriori_1vs1(
      distr_S1 = players[[score[, "joueur_A2"]]][tmp[["keep1"]], ],
      distr_S2 = players[[score[, "joueur_B1"]]][tmp[["keep2"]], ],
      game_len = as.numeric(score[, "game_len"]),
      win = as.numeric(score[, "win"]),
      date = as.Date(score[, "date"]),
      scoreA = as.numeric(scores[, "score_A"]),
      scoreB = as.numeric(scores[, "score_B"]),
      name = c(score[, "joueur_A2"], score[, "joueur_B1"])
    )
    posteriori_per_player <- post_marginal_per_player(posteriori)
    
    posteriori_per_player[[1]][, "p"] <- posteriori_per_player[[1]][, "p"] * (1-probs_ignorees1)
    posteriori_per_player[[2]][, "p"] <- posteriori_per_player[[2]][, "p"] * (1-probs_ignorees2)
    
    players[[score[, "joueur_A2"]]][tmp[["keep1"]], ] <- as.matrix(posteriori_per_player[[1]])
    players[[score[, "joueur_B1"]]][tmp[["keep2"]], ] <- as.matrix(posteriori_per_player[[2]])
    players
  } else {
    
    #très lent donc on va davantage simplifier. discrétisation max de n.
    distrA1 <- distr_simplifier_top_n(players[[score[, "joueur_A1"]]], 10)
    distrA2 <- distr_simplifier_top_n(players[[score[, "joueur_A2"]]], 10)
    distrB1 <- distr_simplifier_top_n(players[[score[, "joueur_B1"]]], 10)
    distrB2 <- distr_simplifier_top_n(players[[score[, "joueur_B2"]]], 10)
    
    players[[score[, "joueur_A1"]]] <- round_up_domain(players[[score[, "joueur_A1"]]])
    players[[score[, "joueur_A2"]]] <- round_up_domain(players[[score[, "joueur_A2"]]])
    players[[score[, "joueur_B1"]]] <- round_up_domain(players[[score[, "joueur_B1"]]])
    players[[score[, "joueur_B2"]]] <- round_up_domain(players[[score[, "joueur_B2"]]])
    
    posteriori <- posteriori_2vs2(
      distr_SA1 = distrA1,
      distr_SA2 = distrA2,
      distr_SB1 = distrB1,
      distr_SB2 = distrB2,
      game_len = as.numeric(score[, "game_len"]),
      win = as.numeric(score[, "win"]),
      date = as.Date(score[, "date"]),
      scoreA = as.numeric(score[, "score_A"]),
      scoreB = as.numeric(score[, "score_B"])
    )
    posteriori_per_player <- post_marginal_per_player(posteriori)
    
    posteriori_per_player[[1]] <- distr_unsimplifier_top_n(distr = posteriori_per_player[[1]], init_distr = players[[score[, "joueur_A1"]]])
    posteriori_per_player[[2]] <- distr_unsimplifier_top_n(distr = posteriori_per_player[[2]], init_distr = players[[score[, "joueur_A2"]]])
    posteriori_per_player[[3]] <- distr_unsimplifier_top_n(distr = posteriori_per_player[[3]], init_distr = players[[score[, "joueur_B1"]]])
    posteriori_per_player[[4]] <- distr_unsimplifier_top_n(distr = posteriori_per_player[[4]], init_distr = players[[score[, "joueur_B2"]]])
    
    players[[score[, "joueur_A1"]]] <- as.matrix(posteriori_per_player[[1]])
    players[[score[, "joueur_A2"]]] <- as.matrix(posteriori_per_player[[2]])
    players[[score[, "joueur_B1"]]] <- as.matrix(posteriori_per_player[[3]])
    players[[score[, "joueur_B2"]]] <- as.matrix(posteriori_per_player[[4]])
    
    players
  }
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
  
  posteriori <- apply(
    matrix(
      apply(scores, 1, function(score) {
      # print(score)
      if(is.na(score["joueur_A1"])) return(likelihood_1vs1_exact(joint_density, score, dataset))
      likelihood_2vs2_exact(joint_density, score, dataset)
    }), nrow = nrow(joint_density$joint_distr), ncol = nrow(scores)
    ),
    1, prod
  ) * joint_density$joint_distr$p
  
  posteriori <- posteriori / sum(posteriori)
  
  #players[pair] <- mapply(simplifier_domain, players[pair], step = ifelse(sapply(players[pair], function(distr) max(distr[, "mu"]) - min(distr[, "mu"])) > 50, 2, 1), SIMPLIFY = FALSE)
  
  #if(max(sapply(players[pair], function(distr) sum(distr[, "p"]))) > 1.0001) stop("Erreur de prob A")

  joint_density$joint_distr$p <- posteriori
  joint_density
}

