source("src/players.R")
source("src/scores.R")
source("src/update_scores.R")
source("src/plots.R")

game_len <- 7


players <- list("a" = matrix(
  c(28, 35, 1, 1, 0.4, 0.6),
  2, dimnames = list(NULL, c("mu", "sig", "p"))
), "b" = matrix(
  c(25, 40, 1, 1, 0.45, 0.55),
  2, dimnames = list(NULL, c("mu", "sig", "p"))
))

distr_mu_sig1 <- players[[1]]
distr_mu_sig2 <- players[[2]]

distr_F1 <- distr_F_finder(distr_mu_sig1)
distr_F2 <- distr_F_finder(distr_mu_sig2)

stopifnot(all.equal(
  sort(apply(distr_F1_F2_1vs1(distr_F1, distr_F2), 1, prod)),
  sort(apply(distr_F1_F2_1vs1(distr_F2, distr_F1), 1, prod))
))

distr_F1_F21 <- distr_F1_F2_1vs1(distr_F1, distr_F2)
distr_F1_F22 <- distr_F1_F2_1vs1(distr_F2, distr_F1)

stopifnot(all.equal(
  sort(distr_P_1vs1(distr_F1_F21)[, "P"]),
  sort(1 - distr_P_1vs1(distr_F1_F22)[, "P"])
))

stopifnot(all.equal(
  {
    distr_P <- distr_P_1vs1(distr_F1_F21)
    distr_P <- cbind(distr_P, "P_win" = p_win_game_of(distr_P[, "P"], game_len))
    sort((distr_P[, "P_win"] * 1 + (1-distr_P[, "P_win"]) * (1-1)) * distr_P[, "p1"] * distr_P[, "p2"])
  }, {
    distr_P <- distr_P_1vs1(distr_F1_F22)
    distr_P <- cbind(distr_P, "P_win" = p_win_game_of(distr_P[, "P"], game_len))
    sort((distr_P[, "P_win"] * 0 + (1-distr_P[, "P_win"]) * (1-0)) * distr_P[, "p1"] * distr_P[, "p2"])
  }
))

stopifnot(all.equal(
  {
    distr_P <- distr_P_1vs1(distr_F1_F21)
    distr_P <- cbind(distr_P, "P_win" = p_win_game_of(distr_P[, "P"], game_len))
    
    Likelihood <- (distr_P[, "P_win"] * 1 + (1-distr_P[, "P_win"]) * (1-1)) * distr_P[, "p1"] * distr_P[, "p2"]
    Likelihood <- Likelihood/sum(Likelihood)
    Likelihood <- cbind(Likelihood, distr_P)
    
    id <- rep(1:(dim_len_mu) * dim_len_F_1vs1, dim_len_mu)
    id <- id + rep((seq(dim_len_mu)-1) * dim_len_mu * dim_len_F_1vs1^2, each = dim_len_mu)
    
    id2 <- rep(rep(rep(1:(dim_len_mu), each=dim_len_F_1vs1), dim_len_F_1vs1), dim_len_mu)
    id2 <- id2 + rep((seq(dim_len_mu)-1) * (dim_len_mu), each = dim_len_mu * dim_len_F_1vs1^2)
    
    posteriori <- cbind(
      Likelihood[id, c("mu1", "sig1", "mu2", "sig2")],
      "p"=sapply(split(Likelihood[, "Likelihood"], id2), sum)
    )
    posteriori[posteriori[, "mu1"] == 35 & posteriori[, "mu2"] == 25, "p"]
  }, {
    sum(Likelihood[Likelihood[, "mu1"] == 35 & Likelihood[, "mu2"] == 25, "Likelihood"])
  }
))



stopifnot(all.equal(
  sort(unlist(posteriori_1vs1(players[[1]], players[[2]], game_len = 7, win = 1))),
  sort(unlist(posteriori_1vs1(players[[2]], players[[1]], game_len = 7, win = 0)))
))



dim_len_F_1vs1 <- 3
dim_len_F_2vs2 <- 3
dim_len_mu <- 2
dim_len_sig <- 1

players <- list("a" = matrix(
  c(28, 35, 1, 1, 0.4, 0.6),
  2, dimnames = list(NULL, c("mu", "sig", "p"))
), "b" = matrix(
  c(25, 40, 1, 1, 0.45, 0.55),
  2, dimnames = list(NULL, c("mu", "sig", "p"))
), "c" = matrix(
  c(16, 20, 1, 1, 0.63, 0.37),
  2, dimnames = list(NULL, c("mu", "sig", "p"))
), "d" = matrix(
  c(51, 49, 1, 1, 0.51, 0.49),
  2, dimnames = list(NULL, c("mu", "sig", "p"))
))

distr_mu_sigA1 <- players[[1]]
distr_mu_sigA2 <- players[[2]]
distr_mu_sigB1 <- players[[3]]
distr_mu_sigB2 <- players[[4]]

distr_FA1 <- distr_F_finder(distr_mu_sigA1, dim_len_F_2vs2)
distr_FA2 <- distr_F_finder(distr_mu_sigA2, dim_len_F_2vs2)
distr_FB1 <- distr_F_finder(distr_mu_sigB1, dim_len_F_2vs2)
distr_FB2 <- distr_F_finder(distr_mu_sigB2, dim_len_F_2vs2)

stopifnot(all.equal(
  sort(apply(distr_F1_F2_2vs2(distr_FA1, distr_FA2, distr_FB1, distr_FB2), 1, sum)),
  sort(apply(distr_F1_F2_2vs2(distr_FB1, distr_FB2, distr_FA2, distr_FA1), 1, sum))
))

distr_F1_F21 <- distr_F1_F2_2vs2(distr_FA1, distr_FA2, distr_FB1, distr_FB2)
distr_F1_F22 <- distr_F1_F2_2vs2(distr_FB1, distr_FB2, distr_FA2, distr_FA1)

stopifnot(all.equal(
  sort(distr_P_2vs2(distr_F1_F21)[, "P"]),
  sort(1 - distr_P_2vs2(distr_F1_F22)[, "P"])
))

stopifnot(all.equal(
  {
    distr_P <- distr_P_2vs2(distr_F1_F21)
    distr_P <- cbind(distr_P, "P_win" = p_win_game_of(distr_P[, "P"], 7))
    sort((distr_P[, "P_win"] * 1 + (1-distr_P[, "P_win"]) * (1-1)) * distr_P[, "pA1"] * distr_P[, "pA2"] * distr_P[, "pB1"] * distr_P[, "pB2"])
  }, {
    distr_P <- distr_P_2vs2(distr_F1_F22)
    distr_P <- cbind(distr_P, "P_win" = p_win_game_of(distr_P[, "P"], 7))
    sort((distr_P[, "P_win"] * 0 + (1-distr_P[, "P_win"]) * (1-0)) * distr_P[, "pA1"] * distr_P[, "pA2"] * distr_P[, "pB1"] * distr_P[, "pB2"])
  }
))

stopifnot(all.equal(
  {
    distr_P <- distr_P_2vs2(distr_F1_F21)
    distr_P <- cbind(distr_P, "P_win" = p_win_game_of(distr_P[, "P"], 7))
    
    Likelihood <- (distr_P[, "P_win"] * 1 + (1-distr_P[, "P_win"]) * (1-1)) * distr_P[, "pA1"] * distr_P[, "pA2"] * distr_P[, "pB1"] * distr_P[, "pB2"]
    Likelihood <- Likelihood/sum(Likelihood)
    Likelihood <- cbind(Likelihood, distr_P)
    
    id <- rep(1:(dim_len_mu) * dim_len_F_2vs2, (dim_len_mu)^3)
    id <- id + rep(rep((seq(dim_len_mu)-1) * dim_len_mu * dim_len_F_2vs2^2, each = dim_len_mu), (dim_len_mu)^2)
    id <- id + rep(rep((seq(dim_len_mu)-1) * (dim_len_mu * dim_len_F_2vs2)^2 * dim_len_F_2vs2, each = (dim_len_mu)^2), dim_len_mu)
    id <- id + rep((seq(dim_len_mu)-1) * (dim_len_mu * dim_len_F_2vs2)^3 * dim_len_F_2vs2, each = (dim_len_mu)^3)
    
    id2 <- rep(rep(1:(dim_len_mu), each=dim_len_F_2vs2), (dim_len_mu * dim_len_F_2vs2)^3)
    id2 <- id2 + rep(rep((seq(dim_len_mu)-1) * (dim_len_mu), each = dim_len_mu * dim_len_F_2vs2^2), (dim_len_mu * dim_len_F_2vs2)^2)
    id2 <- id2 + rep(rep((seq(dim_len_mu)-1) * (dim_len_mu)^2, each = (dim_len_mu * dim_len_F_2vs2)^2 * dim_len_F_2vs2), (dim_len_mu * dim_len_F_2vs2))
    id2 <- id2 + rep((seq(dim_len_mu)-1) * (dim_len_mu)^3, each = (dim_len_mu * dim_len_F_2vs2)^3 * dim_len_F_2vs2)
    
    posteriori <- cbind(
      Likelihood[id, c("muA1", "sigA1", "muA2", "sigA2", "muB1", "sigB1", "muB2", "sigB2")],
      "p"=sapply(split(Likelihood[, "Likelihood"], id2), sum)
    )
    
    posteriori[posteriori[, "muA1"] == 35 & posteriori[, "muA2"] == 25 & posteriori[, "muB1"] == 20 & posteriori[, "muB2"] == 49, "p"]
  }, {
    sum(Likelihood[Likelihood[, "muA1"] == 35 & Likelihood[, "muA2"] == 25 & Likelihood[, "muB1"] == 20 & Likelihood[, "muB2"] == 49, "Likelihood"])
  }
))



stopifnot(all.equal(
  sort(unlist(posteriori_2vs2(players[[1]], players[[2]], players[[3]], players[[4]], game_len = 7, win = 1))),
  sort(unlist(posteriori_2vs2(players[[3]], players[[4]], players[[2]], players[[1]], game_len = 7, win = 0)))
))


dim_len_F_1vs1 <- 1
dim_len_F_2vs2 <- 1
scores <- add_scores(
  data.frame(
    joueur_A1 = c(NA, NA, "d"),
    joueur_A2 = c("a", "a", "b"),
    joueur_B1 = c("b", "b", "c"),
    joueur_B2 = c(NA, NA, "a"),
    win = c(0, 1, 1),
    game_length = c(1, 1, 1)
  ),
  scores_init(),
  date = "2024-03-14"
)

#s'assurer que l'ordre des match n'a pas d'importance
stopifnot(all.equal(
  sapply(update_scores(players, scores[c(1, 2), ]), calculate_ranking),
  sapply(update_scores(players, scores[c(2, 1), ]), calculate_ranking)
))

stopifnot(all.equal(
  sapply(update_scores(players, scores[c(1, 2, 3), ]), calculate_ranking),
  sapply(update_scores(players, scores[c(2, 3, 1), ]), calculate_ranking)
))


scores <- scores_init()
scores <- add_scores(
  data.frame(
    joueur_A1 = c(NA),
    joueur_A2 = c("a", "a"),
    joueur_B1 = c("b", "b"),
    joueur_B2 = c(NA),
    win = c(0, 0),
    game_length = c(7, 6)
  ),
  scores_init(),
  date = "2024-03-14"
)


stopifnot(var(sapply(update_scores(players, scores[1, , drop=F]), calculate_ranking)) > 
            var(sapply(update_scores(players, scores[2, , drop=F]), calculate_ranking)))


                     