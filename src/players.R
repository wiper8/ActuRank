#TODO shiny app?
#TODO se questionner si le score des gens devrait varier si les autres jouent en leur absence.
#TODO vérifier si l'ordonnancement des parties a un impact.


players <- list()

init_distr <- function() {
  mu1 <- qbeta(seq(0, 1, length.out = dim_len_mu), 2, 2) * 99 + 1
  distr_mu1 <- cbind("mu"=mu1, "p"=1/dim_len_mu)
  cbind(mu=mu1, p = 1/dim_len_mu)
}

round_up_domain <- function(distr) {
  distr <- setDT(as.data.frame(distr))
  distr[, mu := round(mu, 1)]
  as.matrix(distr[, .(p = sum(p)), by = mu][order(mu)])
}

tails_simplifier <- function(distr) {
  if(sum(!(distr[, "p"] >= (1 / nrow(distr) / 20) | cumsum(distr[, "p"]) >= 0.001 | cumsum(distr[, "p"]) <= 0.999)) > 0) print("distribution simplifiée")
  distr[, "p"] >= (1 / nrow(distr) / 20) | cumsum(distr[, "p"]) >= 0.001 | cumsum(distr[, "p"]) <= 0.999
}

simplifier_domain <- function(distr, dim_len_mu_min = 15, step = 1) {
  n <- sum(distr[, "p"] >= 1/nrow(distr)/20 | distr[, "p"] >= (max(distr[, "p"]) / 50))
  
  if(n < dim_len_mu_min | sum(head(sort(distr[, "p"], decreasing = T), 15)) > 0.85) {
    distr <- distr_interpolate(distr)
  }
  
  #smooth
  res <- smooth_distr(distr, step = step)
  
  tmp <- tails_simplifier(res)
  
  probs_ignorees <- sum(distr[!tmp, "p"])
  if(probs_ignorees > 0.005) return(res)
  res[tmp, , drop = FALSE]
}

simplifier_joint <- function(joint_density, joint_density_init, seuil = 1 / nrow(joint_density$joint_distr) / 20, max_dimensionality = 20000, absolute_max_dim = 1000000,
                             verbose = FALSE) {
  tmp <- sort(joint_density$joint_distr$p)
  cond_a <- joint_density$joint_distr$p >= seuil
  cond_b <- joint_density$joint_distr$p >= min(tmp[cumsum(tmp) >= 0.005][1], min(tail(tmp, max_dimensionality)))
  cond_c <- joint_density$joint_distr$p >= min(tail(tmp, absolute_max_dim))
  keep <- cond_b & cond_c # cond_a | cond_b
  if(any(!keep) & verbose) {
    print(paste0("% de données conservées : ", round(mean(keep), 5) * 100), collapse = "")
    print(paste0("% de probs conservées : ", round(sum(joint_density$joint_distr$p[keep]), 5) * 100), collapse = "")
  }
  if(sum(joint_density$joint_distr$p[keep]) < 0.95) {
    print(paste0("% de données conservées : ", round(mean(keep), 5) * 100), collapse = "")
    print(paste0("% de probs conservées : ", round(sum(joint_density$joint_distr$p[keep]), 5) * 100), collapse = "")
  }
  
  
  #print(paste0("kept because of seuil : ", round(mean(cond_a[keep]), 3)))
  #print(paste0("kept because of quantile > 0.01 : ", round(mean(cond_b[keep]), 3)))
  #print(paste0("kept because of both : ", round(mean((cond_a | cond_b)[keep]), 3)))
  
  joint_density$joint_distr <- joint_density$joint_distr[keep, ]
  joint_density_init <- joint_density_init[keep]
  joint_density$grid_id <- joint_density$grid_id[keep, ]
  
  joint_density$joint_distr$p <- joint_density$joint_distr$p / sum(joint_density$joint_distr$p)
  joint_density_init <- joint_density_init / sum(joint_density_init)
  list(joint_density, joint_density_init)
}

drift <- function(distr, a = 0.03) {
  #distr * (1-a)^k + priori * (1 - (1-a)^k)
  #\left(1-a\right)^{x}\cdot30\ +\ a\cdot10\cdot\frac{\left(1-\left(1-a\right)^{x}\right)}{a}
  priori <- init_distr()
  x <- c(priori[, "mu"], distr[, "mu"])
  y <- c(a * priori[, "p"], (1-a) * distr[, "p"])
  y <- y[order(x)]
  x <- x[order(x)]
  y <- sapply(unique(x), function(xi) sum(y[x == xi]))
  x <- unique(x)
  cbind(mu=x, p=y)
}


drift_exact <- function(joint_density, joint_density_init, a = 0.03) {
  joint_density$joint_distr$p <- (1 - a) * joint_density$joint_distr$p + a * joint_density_init
  joint_density
}

distr_simplifier_top_n <- function(distr, n = 10) {
  
  if(nrow(distr) <= n) return(distr)
  
  #print("simplifier_top_n")
  
  distr <- rbind(c(distr[1, 1, drop=F], 0), distr) #pour créer une masse à la première valeur
  
  repart <- cumsum(distr[, "p"])
  
  y <- seq(1/(n+1), 1-1/(n+1), length.out=n)
  
  x <- sapply(y, function(p) {
    if(!any(repart <= p)) {
      1 #minimum du range de skill
    } else {
      
    id <- tail(which(repart <= p), 1)
    
    #interpolation linéaire
    (distr[id+1, "mu"] - distr[id, "mu"]) / distr[id + 1, "p"] * (p - repart[id]) + distr[id, "mu"]
    }
  })
  
  y <- rep(1/n, n)
  if(any(x == 1)) {
    y_1 <- sum(y[x == 1])
    y <- y[x != 1]
    x <- x[x != 1]
    x <- c(1, x)
    y <- c(y_1, y)
  }
  res <- matrix(c(x, y), ncol=2, dimnames = dimnames(distr))
  
  #print(ggplot()+
  #  geom_point(aes(x=res[, 1], y=cumsum(res[, 2])))+
  #  geom_line(aes(x=distr[, 1], y=cumsum(distr[, 2])))
  #)
  
  res
}



distr_interpolate <- function(distr) {
  #print("interpolation")
  
  x <- distr[, "mu"]
  y <- distr[, "p"]
  new_x <- (x[-1] + head(x, -1))/2
  new_y <- (y[-1] + head(y, -1))/2
  x <- c(x, new_x)
  y <- c(y, new_y)
  ordre <- order(x)
  y <- y[ordre]
  x <- x[ordre]
  y[c(-1, -length(y))] <- y[c(-1, -length(y))]/sum(y[c(-1, -length(y))]) * (1 - sum(y[c(1, length(y))]))
  res <- cbind(mu=x, p=y)
  smooth_distr(res, step = 1)
}

check_distr <- function(distr) {
  if(nrow(distr) <= 3) warning("domaine très petit")
  if(distr[1, "p"] >= 1/nrow(distr) / 2 & distr[1, "mu"] >= 10) warning("queue gauche plus élevée")
  if(distr[nrow(distr), "p"] >= 1/nrow(distr) / 2 & distr[nrow(distr), "mu"] <= 90) warning("queue droite plus élevée")
}

inverse_cdf <- function(distr) {
  p <- seq(0, 1, 0.001)
  distr <- rbind(c(distr[1, 1] - .Machine$double.eps, p=0), distr)
  distr[, "p"] <- cumsum(distr[, "p"])
  sapply(p, function(p) tail(distr[distr[, "p"] <= p, "mu"], 1))
}

wassertein <- function(distrA, distrB, q = 1) {
  mean((abs(inverse_cdf(distrA) - inverse_cdf(distrB)))^q)^(1 / q)
}

smooth_kernel_distr <- function(distr, init_distr, bandwidth) {
  y <- sapply(init_distr[, "mu"], function(x) sum(distr[, "p"] * dnorm(x, distr[, "mu"], bandwidth)))
  y <- y/sum(y)
  matrix(c(init_distr[, "mu"], y), ncol=2, dimnames = list(NULL, c("mu", "p")))
}

#kernel smoothing
distr_unsimplifier_top_n <- function(distr, init_distr) {
  #print("unsimplifier")
  bandwidth <- optimise(function(bandwidth) {
    kernel <- smooth_kernel_distr(distr, init_distr, bandwidth)
    wassertein(init_distr, kernel)
  }, c(0, 50), tol=0.1)$minimum
  
  smooth_kernel_distr(distr, init_distr, bandwidth)
}

tails_simplifier_1vs1 <- function(distr1, distr2) {
  list(keep1 = tails_simplifier(distr1), keep2 = tails_simplifier(distr2))
}

add_player <- function(name, players) {
  players[[name]] <- init_distr()
  players
}

add_player_exact <- function(name, joint_density, joint_density_init) {
  joint_density$domains[[name]] <- init_distr()[, "mu"]
  
  if(nrow(joint_density$joint_distr) == 0) {
    grid <- as.data.frame(init_distr())
    
    joint_density$grid_id <- data.frame(1:nrow(grid))
    colnames(joint_density$grid_id) <- name
    joint_density_init <- init_distr()[, "p"]
    joint_density$joint_distr <- data.frame(grid["mu"], p = joint_density_init)
    colnames(joint_density$joint_distr) <- c(name, "p")
  } else {
    new_names <- names(joint_density$domains)
    
    tmp <- expand.grid(1:nrow(joint_density$grid_id), init_distr()[, "mu"])
    
    # convert to matrix to speedup
    joint_density$joint_distr <- as.matrix(joint_density$joint_distr)
    joint_density$grid_id <- as.matrix(joint_density$grid_id)
    
    id_col_p <- ncol(joint_density$joint_distr)
    
    grid <- do.call(rbind, lapply(
      1:nrow(tmp), 
      function(i) c(
        joint_density$joint_distr[tmp[i, 1], -id_col_p],
        tmp[i, 2]
      )
    ))
    
    tmp <- expand.grid(1:nrow(joint_density$grid_id), 1:nrow(init_distr()))
    grid_id <- do.call(rbind, lapply(
      1:nrow(tmp), 
      function(i) c(joint_density$grid_id[tmp[i, 1], ], tmp[i, 2])
    ))
    
    joint_density$grid_id <- as.data.frame(grid_id)
    colnames(joint_density$grid_id) <- new_names
    
    joint_density$joint_distr <- as.data.frame(joint_density$joint_distr)
    
    
    joint_density$joint_distr <- cbind(
      as.data.frame(grid),
      p = apply(
        expand.grid(joint_density$joint_distr$p, init_distr()[, "p"]),
        1, prod
      )
    )
    colnames(joint_density$joint_distr) <- c(new_names, "p")
    
    joint_density_init <- apply(
      expand.grid(joint_density_init, init_distr()[, "p"]),
      1, prod
    )
  }
  
  list(joint_density, joint_density_init)
}


compute_credibility <- function(distr, k = 0.1) {
  #e <- sum(distr[, "mu"] * distr[, "p"])
  v <- (sum(distr[, "mu"]^2 * distr[, "p"]) - sum((distr[, "mu"] * distr[, "p"]))^2)
  2 * pnorm(k * 50 / sqrt(v)) - 1
}

is_exact_score_used_for_player <- function(distr, seuil = 0.7) {
  compute_credibility(distr) < seuil
}

