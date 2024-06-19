#TODO ajouter possible drift dans le temps
#TODO shiny app?
#TODO excel partagé et trouver une manière de s'assurer que ce soit de bonnes données
#TODO se questionner si le score des gens devrait varier si les autres jouent en leur absence.
#TODO vérifier si l'ordonnancement des parties a un impact.


players <- list()

dim_len_mu <- 50

init_distr <- function() {
  mu1 <- qbeta(seq(0, 1, length.out = dim_len_mu), 2, 2) * 99 + 1
  distr_mu1 <- cbind("mu"=mu1, "p"=1/dim_len_mu)
  cbind(mu=mu1, p = 1/dim_len_mu)
}

distr_simplifier <- function(distr) {
  #if(sum(!(distr[, "p"] >= (1 / nrow(distr) / 12) | cumsum(distr[, "p"]) >= 0.005 | cumsum(distr[, "p"]) <= 0.995)) > 0) print("distribution simplifiée")
  distr[, "p"] >= (1 / nrow(distr) / 12) | cumsum(distr[, "p"]) >= 0.005 | cumsum(distr[, "p"]) <= 0.995
}

simplifier_domain <- function(distr, dim_len_mu_min = 20) {
  n <- sum(distr[, "p"] >= 1/nrow(distr)/20 | distr[, "p"] >= 1/1000000)
  if(n < dim_len_mu_min | sum(head(sort(distr[, "p"], decreasing = T), 15)) > 0.85) {
    distr <- distr_interpolate(distr)
  }
  #if(sum(!(distr[, "p"] >= 1/nrow(distr)/20 | distr[, "p"] >= 1/1000000)) > 0) print("domain simplified")
  distr[distr[, "p"] >= 1/nrow(distr)/20 | distr[, "p"] >= 1/1000000, ]
}

drift <- function(distr, a = 0.15) {
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
  #print("simplifier_top_n")
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
  print("interpolation")
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

#kernel smoothing
distr_unsimplifier_top_n <- function(distr, init_distr, cap_factor = Inf) {
  #distr <- distr_interpolate(distr_interpolate(distr))
  y <- sapply(init_distr[, "mu"], function(x) sum(distr[, "p"] * dnorm(x, distr[, "mu"], 1.5*mean(diff(distr[, "mu"])))))
  y <- y/sum(y)
  
  factors <- y / init_distr[, "p"]
  #si min est de 0.1 et que j'en ai un autre de 1.2.
  #1.2-1 = 0.2
  factors_change_symmetrical <- ifelse(factors >= 1, factors-1, 1-1/factors)
  factors_change_symmetrical <- pmin(5, abs(factors_change_symmetrical)) * sign(factors_change_symmetrical)
  factors_change_symmetrical <- factors_change_symmetrical * cap_factor / max(abs(factors_change_symmetrical))
  factors <- ifelse(factors_change_symmetrical >= 0, factors_change_symmetrical + 1, 1/(1-factors_change_symmetrical))
  y <- init_distr[, "p"] * factors
  y <- y/sum(y)
  
  #ggplot()+
  #  geom_line(aes(x=0:100, y = sapply(0:100, function(x) sum(dnorm(x, distr[, "mu"], mean(diff(distr[, "mu"])))))))
  
  matrix(c(init_distr[, "mu"], y), ncol=2, dimnames = list(NULL, dimnames(distr)[[2]]))
}

distr_simplifier_1vs1 <- function(distr1, distr2) {
  list(keep1 = distr_simplifier(distr1), keep2 = distr_simplifier(distr2))
}

add_player <- function(name, players) {
  players[[name]] <- init_distr()
  players
}
