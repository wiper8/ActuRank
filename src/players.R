#TODO shiny app?
library(gtools)

players <- list()

init_distr <- function() {
  mu1 <- round(qbeta(seq(0, 1, length.out = dim_len_mu), 1, 1) * 99 + 1, 2)
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

simplifier_joint_dependancy <- function(
    joint_density,
    seuil = 0.0001,
    absolute_max_dim = 1000000,
    min_no_simplif = 5000,
    seuil_lost_lines_pct = 0.05,
    verbose = FALSE,
    count = TRUE
) {
  verbose_overwrite <- verbose
  if (!count) verbose_overwrite <- FALSE
  if (nrow(joint_density$joint_distr) >= min_no_simplif) {
    tmp <- sort(joint_density$joint_distr$p)
    thres_pos <- which(cumsum(tmp) < seuil)
    if(length(thres_pos) != 0) {
      cond_a <- order(order(joint_density$joint_distr$p)) >= max(thres_pos)
    } else cond_a <- TRUE
    cond_b <- joint_density$joint_distr$p >= min(tail(tmp, absolute_max_dim))
    cond_c <- joint_density$joint_distr$p > 0
    cond_d <- order(order(tmp, decreasing = TRUE)) < min_no_simplif
    keep <- (cond_a & cond_b | cond_d) & cond_c
    
    if (mean(keep) >= 1 - seuil_lost_lines_pct) {
      keep <- cond_c
    }
    
    if (count) probs_kept_counter <<- probs_kept_counter * sum(joint_density$joint_distr$p[keep])
    
    if(verbose_overwrite && any(!keep) & verbose) {
      print(paste0("% de données conservées : ", round(mean(keep), 4) * 100), collapse = "")
      print(paste0("% de probs conservées : ", round(sum(joint_density$joint_distr$p[keep]), 6) * 100), collapse = "")
    } else if(verbose_overwrite && sum(joint_density$joint_distr$p[keep]) < 0.99 | sum(keep) > 200000) {
      print(paste0("% de données conservées : ", round(mean(keep), 4) * 100), collapse = "")
      print(paste0("% de probs conservées : ", round(sum(joint_density$joint_distr$p[keep]),6) * 100), collapse = "")
    }
    
  } else {
    keep <- joint_density$joint_distr$p > 0
  }
  
  if(verbose_overwrite && any(!keep) & verbose) {
    print(paste0("nrow joint distr : ", nrow(joint_density$joint_distr[keep, , drop = FALSE]), collapse = ""))
  } else if(verbose_overwrite && sum(joint_density$joint_distr$p[keep]) < 0.99 | sum(keep) > 200000) {
    print(paste0("nrow joint distr : ", nrow(joint_density$joint_distr[keep, , drop = FALSE]), collapse = ""))
  }
  
  joint_density$joint_distr <- joint_density$joint_distr[keep, , drop = FALSE]
  joint_density$grid_id <- joint_density$grid_id[keep, , drop = FALSE]
  joint_density$grid_id <- as.data.frame(
    lapply(as.data.frame(joint_density$grid_id), function(col) as.integer(factor(col)))
  )
  joint_density$joint_distr$p <- joint_density$joint_distr$p / sum(joint_density$joint_distr$p)
  
  joint_density$domains <- lapply(joint_density$joint_distr[1:(ncol(joint_density$joint_distr) - 1)], function(x) sort(unique(x)))
  
  joint_density
}

simplifier_joint_dependancy_for_join_clusters <- function(
    joint_density,
    seuil = 0.0001,
    verbose = TRUE,
    count = TRUE
) {
  verbose_overwrite <- verbose
  if (!count) verbose_overwrite <- FALSE
  tmp <- sort(joint_density$joint_distr$p)
  thres_pos <- which(cumsum(tmp) < seuil)
  if(length(thres_pos) != 0) {
    cond_a <- order(order(joint_density$joint_distr$p)) >= max(thres_pos)
  } else cond_a <- TRUE
  cond_c <- joint_density$joint_distr$p > 0
  keep <- cond_a & cond_c
  
  if (count) probs_kept_counter <<- probs_kept_counter * sum(joint_density$joint_distr$p[keep])
  
  if(verbose_overwrite && any(!keep) & verbose) {
    print(paste0("% de données conservées : ", round(mean(keep), 4) * 100), collapse = "")
    print(paste0("% de probs conservées : ", round(sum(joint_density$joint_distr$p[keep]), 6) * 100), collapse = "")
  } else if(verbose_overwrite && verbose_overwrite && sum(joint_density$joint_distr$p[keep]) < 0.99 | sum(keep) > 200000) {
    print(paste0("% de données conservées : ", round(mean(keep), 4) * 100), collapse = "")
    print(paste0("% de probs conservées : ", round(sum(joint_density$joint_distr$p[keep]),6) * 100), collapse = "")
  }
  
  
  if(any(!keep) & verbose) {
    print(paste0("nrow joint distr : ", nrow(joint_density$joint_distr[keep, , drop = FALSE]), collapse = ""))
  } else if(verbose_overwrite && sum(joint_density$joint_distr$p[keep]) < 0.99 | sum(keep) > 200000) {
    print(paste0("nrow joint distr : ", nrow(joint_density$joint_distr[keep, , drop = FALSE]), collapse = ""))
  }
  
  joint_density$joint_distr <- joint_density$joint_distr[keep, , drop = FALSE]
  joint_density$grid_id <- joint_density$grid_id[keep, , drop = FALSE]
  
  # causerait des erreurs
  # joint_density$grid_id <- as.data.frame(
  #   lapply(as.data.frame(joint_density$grid_id), function(col) as.integer(factor(col)))
  # )
  joint_density$joint_distr$p <- joint_density$joint_distr$p / sum(joint_density$joint_distr$p)
  
  # inutile
  # joint_density$domains <- lapply(joint_density$joint_distr[1:(ncol(joint_density$joint_distr) - 1)], function(x) sort(unique(x)))
  
  joint_density
}

drift_exact <- function(joint_density, clusters, a = 0.05) {
  # trouver la priori marginale
  priori <- new_priori(clusters)
  priori <- priori[priori[, "p"] > 0, , drop = FALSE]
  
  priori_idx <- gtools::combinations(nrow(priori), ncol(joint_density$grid_id), repeats.allowed = TRUE)
  #trouver le domaine plausible de la priori conjointe indépendante du cluster
  priori_p <- apply(priori_idx, 1, function(idx) prod(priori[idx, "p"]))
  # pcq on ne répète pas certaines combins, donc ne somme pas à 1
  # ex 1-2, 2-1, 2-2 (pas un autre 2-2 manquant)
  priori_p_scaled <- priori_p / sum(priori_p)
  
  tmp <- sort(priori_p_scaled)
  thres_pos <- which(cumsum(tmp) < (0.00625 * 2^min(7, ncol(joint_density$grid_id))))
  if(length(thres_pos) != 0) {
    keep <- order(order(priori_p_scaled)) >= max(thres_pos)
  } else {
    keep <- TRUE
  }
  priori_idx <- priori_idx[keep, , drop = FALSE]
  priori_p <- priori_p[keep]
  
  ### ajouter les rows absentes dans joint_density
  # commencer par modifier les grid_id et domaines
  joint_density$grid_id <- as.data.frame(
    mapply(
      function(x) match(x, priori[, "mu"]),
      lapply(seq_len(ncol(joint_density$joint_distr[, -ncol(joint_density$joint_distr), drop = FALSE])), function(i) joint_density$joint_distr[, i]),
      SIMPLIFY = FALSE
    )
  )
  colnames(joint_density$grid_id) <- names(joint_density$domains)
  
  joint_density$domains <- lapply(joint_density$domains, function(x) priori[, "mu"])
  
  # maintenant ajouter les rangées manquantes
  pasted_combins <- apply(joint_density$grid_id, 1, paste0, collapse = "-")
  new_grid_id <- mapply(
    function(idx, new_p) {
      permut <- permute_unique(idx)
      keep <- rep(NA, nrow(permut))
      for (i in 1:ncol(permut)) {
        tmp <- !permut[, i] %in% joint_density$grid_id[, i]
        keep[tmp] <- TRUE
      }
      keep[is.na(keep)] <- apply(
        permut[is.na(keep), , drop = FALSE],
        1,
        function(x) {
          for (i in 1:ncol(permut))
            if (!x[i] %in% joint_density$grid_id[, i]) return(TRUE)
          !any(pasted_combins == paste0(x, collapse = "-"))
        }
      )
      list(permut[keep, , drop = FALSE], new_p)
    }, 
    split(priori_idx, rep(1:nrow(priori_idx), ncol(priori_idx))),
    priori_p,
    SIMPLIFY = FALSE
  )
  
  new_p <- unlist(sapply(new_grid_id, function(ls) rep(ls[[2]], nrow(ls[[1]]))))
  if (length(new_p) > 0) {
    
    
    new_grid_id <- do.call(rbind, sapply(new_grid_id, `[[`, 1))
    new_grid_id <- as.data.frame(new_grid_id)
    colnames(new_grid_id) <- colnames(joint_density$grid_id)
    
    new_distr <- mapply(
      function(i, dom) {
        dom[i]
      },
      lapply(1:ncol(new_grid_id), function(j) new_grid_id[, j]),
      joint_density$domains
    )
    # bugfix
    if (is.null(dim(new_distr))) {
      new_distr <- matrix(new_distr, nrow = 1)
    }
    new_distr <- cbind(
      new_distr,
      p = new_p
    )
    new_distr <- as.data.frame(new_distr)
    colnames(new_distr) <- colnames(joint_density$joint_distr)
    
    
    joint_density$grid_id <- rbind(joint_density$grid_id, new_grid_id)
    joint_density$joint_distr <- rbind(joint_density$joint_distr, cbind(
      new_distr[, -ncol(new_distr), drop = FALSE], p = 0)
    )
  }
  ###
  
  
  priori_likely <- apply(
    joint_density$joint_distr[, -ncol(joint_density$joint_distr), drop = FALSE],
    1,
    function(x) prod(priori[match(x, priori[, "mu"]), "p"])
  )
  
  res <- (1 - a) * joint_density$joint_distr$p + a * priori_likely
  joint_density$joint_distr$p <- res / sum(res)
  joint_density
}

drift_exact2 <- function(joint_density, clusters, a = 0.05) {
  old_dom <- joint_density$domains
  
  # trouver la priori marginale pour chaque joueur centrée a l'espérance
  priori <- new_priori_centered(joint_density)
  
  priori_idx_nrow <- prod(sapply(priori, function(x) nrow(x)))
  
  #trouver les probs de la priori conjointe indépendante du cluster
  if (!all(sapply(priori, function(distr) length(unique(distr[, "p"])) == 1))) stop("drift non uniforme non implémenté")
  priori_p <- rep(prod(sapply(priori, function(x) x[1, "p"])), priori_idx_nrow)
  
  joint_density$domains <- mapply(
    function(dom, prio) sort(unique(c(prio[, "mu"], dom))),
    joint_density$domains, 
    priori, 
    SIMPLIFY = FALSE
  )
  
  ### ajouter les rows absentes dans joint_density
  # commencer par modifier les grid_id et domaines
  joint_density$grid_id <- as.data.frame(
    mapply(
      function(x, dom) match(x, dom),
      lapply(seq_len(ncol(joint_density$joint_distr[, -ncol(joint_density$joint_distr), drop = FALSE])), function(i) joint_density$joint_distr[, i]),
      joint_density$domains,
      SIMPLIFY = FALSE
    )
  )
  colnames(joint_density$grid_id) <- names(joint_density$domains)
  
  # maintenant ajouter les rangées manquantes
  pasted_combins <- apply(joint_density$grid_id, 1, paste0, collapse = "-")
  priori_idx <- expand.grid(lapply(seq_len(ncol(joint_density$grid_id)), function(j) {
    match(priori[[j]][, "mu"], joint_density$domains[[j]])
  }))
  colnames(priori_idx) <- names(joint_density$domains)
  pasted_combin_new <- apply(priori_idx, 1, paste0, collapse = "-")
  
  keep <- !pasted_combin_new %in% pasted_combins
  priori_idx <- priori_idx[keep, , drop = FALSE]
  priori_p_unif <- priori_p[1]
  priori_p <- priori_p[keep]
  
  if (length(priori_p) > 0) {
    new_distr <- mapply(
      function(i, dom) {
        dom[i]
      },
      lapply(1:ncol(priori_idx), function(j) priori_idx[, j]),
      joint_density$domains
    )
    # bugfix
    if (is.null(dim(new_distr))) {
      new_distr <- matrix(new_distr, nrow = 1)
    }
    new_distr <- cbind(
      new_distr,
      p = priori_p
    )
    new_distr <- as.data.frame(new_distr)
    colnames(new_distr) <- colnames(joint_density$joint_distr)
    
    joint_density$grid_id <- rbind(joint_density$grid_id, priori_idx)
    joint_density$joint_distr <- rbind(joint_density$joint_distr, cbind(
      new_distr[, -ncol(new_distr), drop = FALSE], p = 0)
    )
  }
  
  # identifier les combinaisons qui ne sont pas touchées par le drift
  old_idx_not_in_drift_priori_joint <- mapply(
    function(prio, dom, new_dom) {
      #tmp <- match(prio[, "mu"], dom)
      #seq_along(dom)[!seq_along(dom) %in% tmp]
      
      id <- seq_along(new_dom)
      tmp <- setdiff(dom, prio[, "mu"])
      id[new_dom %in% tmp]
    },
    priori, old_dom, joint_density$domains,
    SIMPLIFY = FALSE
  )
  touched_by_new_drift_joint <- !apply(
    mapply(function(id, id_new) id %in% id_new, joint_density$grid_id, old_idx_not_in_drift_priori_joint),
    1, any
  )
  
  res <- (1 - a) * joint_density$joint_distr$p + a * priori_p_unif * touched_by_new_drift_joint
  joint_density$joint_distr$p <- res / sum(res)
  joint_density
}

drift_exact3 <- function(joint_density, clusters, a = 0.01) {
  joint_density$joint_distr$p <- (1 - a) * joint_density$joint_distr$p + a / nrow(joint_density$joint_distr)
  joint_density
}

test_hyp <- function(clusters) {
  # sorter les joueurs par classement
  players <- marginal_from_joint_dependancy(clusters)
  ranks <- sapply(players, function(distr) calculate_skill(distr, players))
  ordre <- order(ranks, decreasing = TRUE)
  players <- players[ordre]
  ranks <- ranks[ordre]
  # faire la grille de tous les paires possibles
  pairs <- t(combn(1:length(ranks), 2))
  
  res <- matrix(
    c(
      names(players)[pairs[, 1]],
      names(players)[pairs[, 2]],
      rep(NA, nrow(pairs)),
      rep(NA, nrow(pairs))
    ), ncol = 4,
    dimnames = list(NULL, c("A", "B", "P(A>B)", "P(A=B)"))
  )
  
  # pour chaque paire :
  for (i in 1:nrow(pairs)) {
    groups_to_join <- which(sapply(lapply(clusters, `[[`, "names"), function(noms) any(names(players)[pairs[i, ]] %in% noms)))
    # si pas dans le meme clusters, 
    if (length(groups_to_join) > 1) {
      # expand.grid les 2 marginales en X les probs
      A <- players[[pairs[i, 1]]]
      B <- players[[pairs[i, 2]]]
      joint <- data.frame(rep(A[, "mu"], nrow(B)), rep(B[, "mu"], each = nrow(A)), p = rep(A[, "p"], nrow(B)) * rep(B[, "p"], each = nrow(A)))
      colnames(joint) <- c(names(players)[pairs[i, ]], "p")
    } else {
      joint <- clusters[[groups_to_join]]$joint_distr
    }
    # sommes la joint_distr par paire de joueurs
    subset <- joint[, c(names(players)[pairs[i, ]], "p")]
    res[i, 3] <- sum(apply(subset, 1, function(x) (x[1] > x[2]) * x[3]))
    res[i, 4] <- sum(apply(subset, 1, function(x) (x[1] == x[2]) * x[3]))
  }
  
  res <- as.data.frame(res)
  res$`P(A>B)` <- as.numeric(res$`P(A>B)`)
  res$`P(A=B)` <- as.numeric(res$`P(A=B)`)
  res
}

# distr_simplifier_top_n <- function(distr, n = 10) {
#   
#   if(nrow(distr) <= n) return(distr)
#   
#   #print("simplifier_top_n")
#   
#   distr <- rbind(c(distr[1, 1, drop=F], 0), distr) #pour créer une masse à la première valeur
#   
#   repart <- cumsum(distr[, "p"])
#   
#   y <- seq(1/(n+1), 1-1/(n+1), length.out=n)
#   
#   x <- sapply(y, function(p) {
#     if(!any(repart <= p)) {
#       1 #minimum du range de skill
#     } else {
#       
#     id <- tail(which(repart <= p), 1)
#     
#     #interpolation linéaire
#     (distr[id+1, "mu"] - distr[id, "mu"]) / distr[id + 1, "p"] * (p - repart[id]) + distr[id, "mu"]
#     }
#   })
#   
#   y <- rep(1/n, n)
#   if(any(x == 1)) {
#     y_1 <- sum(y[x == 1])
#     y <- y[x != 1]
#     x <- x[x != 1]
#     x <- c(1, x)
#     y <- c(y_1, y)
#   }
#   res <- matrix(c(x, y), ncol=2, dimnames = dimnames(distr))
#   
#   #print(ggplot()+
#   #  geom_point(aes(x=res[, 1], y=cumsum(res[, 2])))+
#   #  geom_line(aes(x=distr[, 1], y=cumsum(distr[, 2])))
#   #)
#   
#   res
# }


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

# smooth_kernel_distr <- function(distr, init_distr, bandwidth) {
#   y <- sapply(init_distr[, "mu"], function(x) sum(distr[, "p"] * dnorm(x, distr[, "mu"], bandwidth)))
#   y <- y/sum(y)
#   matrix(c(init_distr[, "mu"], y), ncol=2, dimnames = list(NULL, c("mu", "p")))
# }

#kernel smoothing
# distr_unsimplifier_top_n <- function(distr, init_distr) {
#   #print("unsimplifier")
#   bandwidth <- optimise(function(bandwidth) {
#     kernel <- smooth_kernel_distr(distr, init_distr, bandwidth)
#     wassertein(init_distr, kernel)
#   }, c(0, 50), tol=0.1)$minimum
#   
#   smooth_kernel_distr(distr, init_distr, bandwidth)
# }

# tails_simplifier_1vs1 <- function(distr1, distr2) {
#   list(keep1 = tails_simplifier(distr1), keep2 = tails_simplifier(distr2))
# }

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

add_player_dependancy <- function(name, clusters) {
  if (length(clusters) == 0) {
    grid <- as.data.frame(init_distr())
    
    grid_id <- data.frame(1:nrow(grid))
    colnames(grid_id) <- name
    joint_density_init <- init_distr()[, "p"]
    joint_distr <- data.frame(grid["mu"], p = joint_density_init)
    colnames(joint_distr) <- c(name, "p")
    
    dom <- list(init_distr()[, "mu"])
    names(dom) <- name
    
    clusters[[length(clusters) + 1]] <- list(
      grid_id = grid_id,
      joint_distr = joint_distr,
      domains = dom,
      names = name
    )
    return(clusters)
  }
  priori <- new_priori(clusters)
  priori <- priori[priori[, "p"] > 0, ]
  grid <- as.data.frame(priori)
  
  grid_id <- data.frame(1:nrow(grid))
  colnames(grid_id) <- name
  joint_density_init <- priori[, "p"]
  joint_distr <- data.frame(grid["mu"], p = joint_density_init)
  colnames(joint_distr) <- c(name, "p")
  
  dom <- list(priori[, "mu"])
  names(dom) <- name
  
  clusters[[length(clusters) + 1]] <- list(
    grid_id = grid_id,
    joint_distr = joint_distr,
    domains = dom,
    names = name
  )
  clusters
}

new_priori <- function(clusters) {
  
  marginales <- marginal_from_joint_dependancy(clusters)
  x <- init_distr()[, "mu"]
  matrix(
    c(x,
      sapply(x, function(x_i) {
        mean(sapply(marginales, function(distr) {
          if (any(distr[, "mu"] == x_i)) {
            distr[distr[, "mu"] == x_i, "p"]
          } else {
            0
          }
        }))
      })),
    ncol = 2,
    dimnames = dimnames(init_distr())
  )
}

new_priori_centered <- function(joint_density) {
  
  marginales <- marginal_from_joint_dependancy(list(joint_density))
  priori <- init_distr()
  if (length(unique(priori[, "p"])) != 1) stop("non implémenté pour priori non uniforme")
  bounds <- range(sapply(marginales, function(x) range(x[, "mu"])))
  priori <- priori[priori[, "mu"] >= bounds[1] & priori[, "mu"] <= bounds[2], , drop = FALSE]
  lapply(marginales, function(x) {
    mean_mu <- sum(x[, "mu"] * x[, "p"])
    range_to_closest_bound <- min(c(abs(bounds[1] - mean_mu), abs(bounds[2] - mean_mu)))
    new_domain <- priori[mean_mu - range_to_closest_bound <= priori[, "mu"] & mean_mu + range_to_closest_bound >= priori[, "mu"], "mu"]
    matrix(
      c(
        new_domain,
        rep(1 / length(new_domain), length(new_domain))
      ),
      ncol = 2,
      dimnames = list(NULL, c("mu", "p"))
    )  
  })
}

join_clusters <- function(clusters, which_i, seuil_simplif = 0.00001, count = TRUE) {
  if (length(which_i) == 1) {
    return(clusters[[which_i]])
  }
  
  clusters_subset <- clusters[which_i]
  
  new_names <- unlist(lapply(clusters_subset, `[[`, "names"))
  
  new_grid_id <- do.call(
    expand.grid,
    lapply(clusters_subset, function(clust)
      1:nrow(clust$grid_id)
    )
  )
  
  if (nrow(new_grid_id) > 100000) print(paste0("joining big grid > 100k (", nrow(new_grid_id), ")"))
  
  id_col_p <- sapply(
    clusters_subset,
    function(clust) ncol(clust$joint_distr)
  )
  
  # convert to matrix to speedup
  new_grid_id <- as.matrix(new_grid_id)
  clusters_subset <- lapply(clusters_subset, function(clust) {
    clust$joint_distr <- as.matrix(clust$joint_distr)
    clust
  })
  
  prob_joint <- sapply(
    1:nrow(new_grid_id),
    function(i) {
      a <- new_grid_id[i, ]
      prod(mapply(function(clust, j, id_col) clust$joint_distr[j, id_col], clusters_subset, a, id_col_p))
    }
  )
  
  # simplifier un peu les cas quasi-impossible pour grandement accélérer
  tmp <- simplifier_joint_dependancy_for_join_clusters(
    joint_density = list(
      grid_id = new_grid_id,
      joint_distr = data.frame(A = rep(NA, length(prob_joint)), p = prob_joint)
    ),
    seuil = seuil_simplif,
    verbose = TRUE,
    count = count
  )
  
  new_grid_id <- tmp$grid_id
  prob_joint <- tmp$joint_distr$p
  ###
  
  
  grid <- do.call(
    rbind,
    lapply(
      1:nrow(new_grid_id),
      function(i) {
        a <- new_grid_id[i, ]
        unlist(mapply(function(clust, j, id_col) clust$joint_distr[j, -id_col, drop = FALSE], clusters_subset, a, id_col_p, SIMPLIFY = FALSE))
      }
    )
  )
  
  new_grid_id_wide <- do.call(
    rbind,
    lapply(
      1:nrow(new_grid_id),
      function(i) {
        a <- new_grid_id[i, ]
        unlist(mapply(function(clust, j) clust$grid[j, , drop = FALSE], clusters_subset, a, SIMPLIFY = FALSE))
      }
    )
  )
  
  new_joint_distr <- as.data.frame(cbind(grid, p = prob_joint))
  
  new_grid_id_wide <- as.data.frame(new_grid_id_wide)
  
  colnames(new_grid_id_wide) <- new_names
  colnames(new_joint_distr) <- c(new_names, "p")
  
  list(
    grid_id = new_grid_id_wide,
    joint_distr = new_joint_distr,
    domains = unlist(lapply(clusters_subset, function(clust) clust$domains), recursive = FALSE)
  )
}

join_clusters_pair_marginal <- function(clusters, names) {
  which_i <- sapply(lapply(clusters, `[[`, "names"), function(noms) any(noms %in% names))
  
  # summarise marginal
  marginals <- lapply(clusters[which_i], function(clust) {
    if (length(clust$names) == 1) return(clust)
    players_to_marginalize <- names[which(names %in% clust$names)]
    which_out <- which(colnames(clust$grid_id) %in% players_to_marginalize)
    if (isTRUE(all.equal(seq_along(clust$grid_id), sort(which_out)))) {
      clust
    } else {
      marginal_joint_dependancy_out_only(clust, which_out)
    }
  })
  
  join_clusters(marginals, seq_along(marginals), seuil_simplif = 0, count = FALSE)
}

generate_partitions <- function(n, k, min_cluster_size = 1) {
  all_partitions <- lapply(1:ceiling((n - 1) / 2), function(k) {
    if (n %% 2 == 0 & n / 2 == k) {
      tmp <- combn(1:n, k, simplify = FALSE)
      tmp[sapply(tmp, function(x) 1 %in% x)]
    } else {
      combn(1:n, k, simplify = FALSE)
    }
  })
  
  # partitions <- lapply(all_partitions, function(combs) {
  #   lapply(combs, function(cluster1) {
  #     cluster2 <- setdiff(1:n, cluster1)
  #     list(cluster1, cluster2)
  #   })
  # })
  
  partitions <- unlist(all_partitions, recursive = FALSE)
  
  partitions <- lapply(partitions, function(x) if(min(length(x), (n - length(x))) >= min_cluster_size) x else NULL)
  partitions <- partitions[!sapply(partitions, is.null)]
  partitions
}

permutations <- function(n, v) {
  sub <- function(n, v) {
    if (n == 1) 
      matrix(v, 1, 1)
    else {
      X <- NULL
      for (i in 1:n) X <- rbind(X, cbind(v[i], Recall(n - 1, v[-i])))
      X
    }
  }
  sub(n, v)
}

permute_unique <- function(nums) {
  # Sort the input vector to ensure duplicates are adjacent
  nums <- sort(nums)
  
  # List to store the final result
  result <- list()
  
  # Backtracking function to build permutations
  backtrack <- function(path, remaining) {
    if (length(remaining) == 0) {
      # Append current permutation to the result list
      result[[length(result) + 1]] <<- path
      return()
    }
    
    for (i in seq_along(remaining)) {
      # Skip duplicates by checking if the current element is the same as the previous one
      if (i > 1 && remaining[i] == remaining[i - 1]) {
        next
      }
      
      # Recursively backtrack with the current element added to the path
      backtrack(c(path, remaining[i]), remaining[-i])
    }
  }
  
  # Start the backtracking process with an empty path and the sorted input vector
  backtrack(c(), nums)
  
  do.call(rbind, result)
}


recluster_dependancy <- function(joint_density, dataset, excluded_members = NULL,
                                 max_cluster_size = Inf, min_cluster_size = 1,
                                 joint_distr_size_skip = 500) {
  if (dataset == "ping") MI_thresh_for_indep <- 0.01
  else MI_thresh_for_indep <- 0
  
  if (
    length(joint_density$domains) == 1 ||
    nrow(joint_density$joint_distr) <= joint_distr_size_skip ||
    MI_thresh_for_indep == 0
  ) {
    joint_density$names <- names(joint_density$domains)
    return(list(list(joint_density), excluded_members))
  }
  
  # faire un gros cluster et pour chq élément regarder si on peut les retirer d'une manière greedy
  combins <- generate_partitions(length(joint_density$domains), 2, min_cluster_size = min_cluster_size)
  combins_len <- sapply(combins, length)
  
  if (length(combins) > 0) {
    MI <- rep(Inf, length(combins))
    exclude_combins_improbable <- c()
    # faire toutes les combins qu'on rejette 1 personne, ensuite toutes rejet 2pers, etc.
    for (len_out in unique(combins_len)) {
      for (j in which(combins_len == len_out)) {
        if (!any(combins[[j]] %in% exclude_combins_improbable) & !any(joint_density$names[combins[[j]]] %in% excluded_members)) {
          MI[j] <- relative_mutual_information(joint_density, combins[[j]])
          if (MI[j] > 0.8 * len_out / length(joint_density$domains)) {
            print(paste0("MI : ", round(MI[j], 3)))
            excluded_members <- unique(c(excluded_members, joint_density$names[combins[[j]]]))
            print(paste0("excluded members : ", paste0(excluded_members, ", "), collapse = ""))
            exclude_combins_improbable <- unique(c(combins[[j]], exclude_combins_improbable))
          }
          if (MI[j] > 0.5 * len_out / length(joint_density$domains)) {
            exclude_combins_improbable <- unique(c(combins[[j]], exclude_combins_improbable))
          }
        }
      }
      
      # si le min de MI de ce groupe est sous le seuil, on essaie d'augmenter la 
      # taille de rejetés (donc on continue à chercher plus de combins)
      # pcq plus on rejette bcp, plus min(MI) sera élevé donc on veut arrêter
      # quand tout le groupe dépasse le seuil
      if (
        min(MI[which(combins_len == len_out)]) >= MI_thresh_for_indep ||
        (len_out > 1 && (
          min(MI[which(combins_len == len_out - 1)]) < min(MI[which(combins_len == len_out)])
        ))
      ){
        # print(paste0("MI : ", paste0(round(MI, 3), collapse = ", ")))
        # on compare toujours un séparation entre 2 subsets, comme si on a 2 variables
        # MI = H(X) + H(Y) - H(X,Y)
        # max(H(X), H(Y)) <= H(X)+H(Y)-max(H(X), H(Y))
        # borne sup(MI) = H(X) + H(Y) - max(H(x),H(Y))
        if (min(MI) < MI_thresh_for_indep) {
          print(paste0(sum(MI < Inf), " / ", length(MI)))
          print(paste0("MI : ", round(min(MI), 5), collapse = ""))
          # l'exclure du gros cluster in
          # TODO moyen d'optimiser car dédouble les calculs déjà faits dans mutual_info()
          new_clusters <- marginal_joint_dependancy(joint_density, combins[[which.min(MI)]], format = 2)
          tmp <- recluster_dependancy(new_clusters[[1]], dataset, excluded_members)
          return(
            list(
              c(
                new_clusters[2],
                tmp[[1]]
              ),
              tmp[[2]]
            )
          )
        } else {
          break
        }
      }
    }
    if (length(joint_density$domains) > max_cluster_size) {
      MI <- MI[length(joint_density$domains) - sapply(combins, length) <= max_cluster_size]
      combins <- combins[length(joint_density$domains) - sapply(combins, length) <= max_cluster_size]
      new_clusters <- marginal_joint_dependancy(
        joint_density,
        combins[[which.min(MI)]],
        format = 2
      )
      tmp <- recluster_dependancy(
        new_clusters[[1]],
        dataset,
        excluded_members
      )
      return(
        list(
          c(
            new_clusters[2],
            tmp[[1]]
          ),
          tmp[[2]]
        )
      )
    }
  }
  joint_density$names <- unlist(names(joint_density$domains))
  list(list(joint_density), excluded_members)
}

marginal_joint_dependancy <- function(joint_density, idx_out, format = 1) {
  new_grid_in <- as.matrix(joint_density$grid_id[, -idx_out, drop = FALSE])
  new_grid_out <- as.matrix(joint_density$grid_id[, idx_out, drop = FALSE])
  
  new_values_grid_in <- unique(lapply(1:nrow(new_grid_in), function(i) new_grid_in[i, , drop = TRUE]))
  new_values_grid_out <- unique(lapply(1:nrow(new_grid_out), function(i) new_grid_out[i, , drop = TRUE]))
  
  marginal_grid_in <- as.data.frame(do.call(
    rbind,
    new_values_grid_in
  ))
  marginal_grid_out <- as.data.frame(do.call(
    rbind,
    new_values_grid_out
  ))
  
  a_in <- sapply(new_values_grid_in, paste, collapse = "-")
  b_in <- apply(new_grid_in, 1, paste, collapse = "-")
  a_out <- sapply(new_values_grid_out, paste, collapse = "-")
  b_out <- apply(new_grid_out, 1, paste, collapse = "-")
  
  b_in <- factor(b_in, levels = a_in)
  b_out <- factor(b_out, levels = a_out)
  
  marginal_prob_in <- as.vector(tapply(joint_density$joint_distr$p, b_in, sum))
  marginal_prob_out <- as.vector(tapply(joint_density$joint_distr$p, b_out, sum))
  
  if (format == 1) {
    return(list(
      marginal_prob_in = marginal_prob_in,
      new_values_grid_in = new_values_grid_in,
      marginal_prob_out = marginal_prob_out,
      new_values_grid_out = new_values_grid_out
    ))
  }
  
  distr_in <- do.call(rbind,
                      lapply(1:nrow(marginal_grid_in), function(i) {
                        sapply(1:ncol(marginal_grid_in), function(j) {
                          joint_density$domains[-idx_out][[j]][marginal_grid_in[i, j]]
                        })
                      })
  )
  distr_out <- do.call(rbind,
                       lapply(1:nrow(marginal_grid_out), function(i) {
                         sapply(1:ncol(marginal_grid_out), function(j) {
                           joint_density$domains[idx_out][[j]][marginal_grid_out[i, j]]
                         })
                       })
  )
  distr_in <- as.data.frame(distr_in)
  distr_in <- cbind(distr_in, marginal_prob_in)
  colnames(distr_in) <- c(names(joint_density$domains)[-idx_out], "p")
  distr_out <- as.data.frame(distr_out)
  distr_out <- cbind(distr_out, marginal_prob_out)
  colnames(distr_out) <- c(names(joint_density$domains)[idx_out], "p")
  
  list(
    list(
      grid_id = marginal_grid_in,
      joint_distr = distr_in,
      domains = lapply(distr_in[1:(ncol(distr_in) - 1)], function(x) sort(unique(x))),
      names = colnames(distr_in)[-ncol(distr_in)]
    ),
    list(
      grid_id = marginal_grid_out,
      joint_distr = distr_out,
      domains = lapply(distr_out[1:(ncol(distr_out) - 1)], function(x) sort(unique(x))),
      names = colnames(distr_out)[-ncol(distr_out)]
    )
  )
}

marginal_joint_dependancy_out_only <- function(joint_density, idx_out) {
  new_grid_out <- as.matrix(joint_density$grid_id[, idx_out, drop = FALSE])
  
  new_values_grid_out <- unique(lapply(1:nrow(new_grid_out), function(i) new_grid_out[i, , drop = TRUE]))
  
  marginal_grid_out <- as.data.frame(do.call(
    rbind,
    new_values_grid_out
  ))
  
  a_out <- sapply(new_values_grid_out, paste, collapse = "-")
  b_out <- apply(new_grid_out, 1, paste, collapse = "-")
  
  b_out <- factor(b_out, levels = a_out)
  
  marginal_prob_out <- as.vector(tapply(joint_density$joint_distr$p, b_out, sum))
  
  distr_out <- do.call(rbind,
                       lapply(1:nrow(marginal_grid_out), function(i) {
                         sapply(1:ncol(marginal_grid_out), function(j) {
                           joint_density$domains[idx_out][[j]][marginal_grid_out[i, j]]
                         })
                       })
  )
  distr_out <- as.data.frame(distr_out)
  distr_out <- cbind(distr_out, marginal_prob_out)
  colnames(distr_out) <- c(names(joint_density$domains)[idx_out], "p")
  
  list(
      grid_id = marginal_grid_out,
      joint_distr = distr_out,
      domains = lapply(distr_out[1:(ncol(distr_out) - 1)], function(x) sort(unique(x))),
      names = colnames(distr_out)[-ncol(distr_out)]
    )
}

entropy <- function(p) -sum(p * log(p, 2))

relative_mutual_information <- function(joint_density, idx_out) {
  tmp <- marginal_joint_dependancy(joint_density, idx_out, format = 1)
  marginal_prob_in = tmp$marginal_prob_in
  new_values_grid_in = tmp$new_values_grid_in
  marginal_prob_out = tmp$marginal_prob_out
  new_values_grid_out = tmp$new_values_grid_out
  
  new_values_grid_in <- do.call(rbind, new_values_grid_in)
  new_values_grid_out <- do.call(rbind, new_values_grid_out)
  
  
  a_in <- apply(joint_density$grid_id[, -idx_out, drop = FALSE], 1, paste, collapse = "-")
  b_in <- apply(new_values_grid_in, 1, paste, collapse = "-")
  a_out <- apply(joint_density$grid_id[, idx_out, drop = FALSE], 1, paste, collapse = "-")
  b_out <- apply(new_values_grid_out, 1, paste, collapse = "-")
  
  deno_in <- marginal_prob_in[match(a_in, b_in)]
  deno_out <- marginal_prob_out[match(a_out, b_out)]
  
  entropy_in <- entropy(marginal_prob_in)
  entropy_out <- entropy(marginal_prob_in)
  
  mutual_info <- sum(joint_density$joint_distr$p * (log(joint_density$joint_distr$p, 2) - log(deno_in, 2) - log(deno_out, 2)))
  res <- mutual_info / (entropy_in + entropy_out - max(entropy_in, entropy_out))
  if (mutual_info == 0 | mutual_info < 10^-10 | is.nan(mutual_info)) {
    res <- 1
  }
  if (!is.numeric(res) | is.na(res)) {
    print(mutual_info)
    print(entropy_in)
    print(entropy_out)
  }
  if (res > 1.00001) {
    print(res)
    print(mutual_info)
    print(entropy_in)
    print(entropy_out)
    stop("Mutual information / max  > 1 shouldn't be possible")
  }
  res
}


compute_credibility <- function(distr, k = 100 / dim_len_mu / 50) {
  # crédibilité plus traditionnelle
  # e <- sum(distr[, "mu"] * distr[, "p"])
  # v <- (sum(distr[, "mu"]^2 * distr[, "p"]) - e^2)
  # v <- max(v, 0) # car imprécisions sur les opérations floating point
  # 2 * pnorm(k * 50 / sqrt(v)) - 1 #2 * pnorm(k * e / sqrt(v)) - 1
  
  # crédibilitée basée sur la shape dla distribution
  # P(-50k <= X-E[X] <= 50k)
  e <- sum(distr[, "mu"] * distr[, "p"])
  sum(distr[(distr[, "mu"] - e) <= 50*k & -50*k <= (distr[, "mu"] - e), "p"])
}

compute_multivariate_credibility <- function(distr) {
  p_ncol <- ncol(distr)
  idx <- t(apply(distr[, -p_ncol, drop = FALSE], 1, order2))
  spread <- diff(range(idx))
  apply(idx, 2, function(id_row) {
    
    e <- sum(id_row * distr[, p_ncol])
    e2 <- sum(id_row^2 * distr[, p_ncol])
    v <- e2 - e^2
    # v=0 -> credib = 1
    # v = ((b-a+1)^2 - 1)/12 -> credib = 0
    1 + (-12/((spread + 1)^2 - 1)) * v
  })
}

order2 <- function(x) {
  unique_values <- sort(unique(x))
  sapply(x, function(x) which(x == unique_values))
}

is_exact_score_used_for_player <- function(distr, seuil = 0.95) {
  compute_multivariate_credibility(distr) < seuil
}

