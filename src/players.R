#TODO shiny app?
#TODO se questionner si le score des gens devrait varier si les autres jouent en leur absence.
#TODO vérifier si l'ordonnancement des parties a un impact.


players <- list()

dim_len_mu <- 30

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

distr_simplifier <- function(distr) {
  if(sum(!(distr[, "p"] >= (1 / nrow(distr) / 12) | cumsum(distr[, "p"]) >= 0.005 | cumsum(distr[, "p"]) <= 0.995)) > 0) print("distribution simplifiée")
  distr[, "p"] >= (1 / nrow(distr) / 12) | cumsum(distr[, "p"]) >= 0.005 | cumsum(distr[, "p"]) <= 0.995
}

simplifier_domain <- function(distr, dim_len_mu_min = 20) {
  n <- sum(distr[, "p"] >= 1/nrow(distr)/20 | distr[, "p"] >= (max(distr[, "mu"]) / 50))
  
  if(n < dim_len_mu_min | sum(head(sort(distr[, "p"], decreasing = T), 15)) > 0.85) {
    distr <- distr_interpolate(distr)
  }
  
  #remove too low points
  distr <- distr[distr[, "p"] >= 1/nrow(distr)/20 | distr[, "p"] >= (max(distr[, "mu"]) / 50), , drop = FALSE]
  
  #smooth
  smooth_distr(distr, step = 2)
  
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
  cbind(mu=x, p=y)
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
  cbind(mu=init_distr[, "mu"], p=y)
}

#kernel smoothing
distr_unsimplifier_top_n <- function(distr, init_distr) {
  #print("unsimplifier")
  bandwidth <- optimise(function(bandwidth) {
    kernel <- smooth_kernel_distr(distr, init_distr, bandwidth)
    wassertein(init_distr, kernel)
    
  }, c(0, 100), tol=0.1)$minimum
  
  smooth_kernel_distr(distr, init_distr, bandwidth)
}

distr_simplifier_1vs1 <- function(distr1, distr2) {
  list(keep1 = distr_simplifier(distr1), keep2 = distr_simplifier(distr2))
}

add_player <- function(name, players) {
  players[[name]] <- init_distr()
  players
}


compute_credibility <- function(distr, k = 0.1) {
  #e <- sum(distr[, "mu"] * distr[, "p"])
  v <- (sum(distr[, "mu"]^2 * distr[, "p"]) - sum((distr[, "mu"] * distr[, "p"]))^2)
  2 * pnorm(k * 50 / sqrt(v)) - 1
}

is_exact_score_used_for_player <- function(distr, seuil = 0.7) {
  compute_credibility(distr) < seuil
}

