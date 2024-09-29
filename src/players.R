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

simplifier_joint_dependancy <- function(
    joint_density,
    seuil = 0.0001,
    absolute_max_dim = 1000000,
    min_no_simplif = 5000,
    verbose = FALSE
) {
  if (nrow(joint_density$joint_distr) >= min_no_simplif) {
    tmp <- sort(joint_density$joint_distr$p)
    thres_pos <- max(which(cumsum(tmp) < seuil))
    if(length(thres_pos) != 0) {
      cond_a <- order(order(joint_density$joint_distr$p)) >= thres_pos # min(thres_pos, length(joint_density$joint_distr$p) - min_no_simplif)
    } else cond_a <- TRUE
    cond_b <- joint_density$joint_distr$p >= min(tail(tmp, absolute_max_dim))
    cond_c <- joint_density$joint_distr$p > 0
    keep <- cond_a & cond_b & cond_c
    
    if(any(!keep) & verbose) {
      print(paste0("% de données conservées : ", round(mean(keep), 4) * 100), collapse = "")
      print(paste0("% de probs conservées : ", round(sum(joint_density$joint_distr$p[keep]), 4) * 100), collapse = "")
    } else if(sum(joint_density$joint_distr$p[keep]) < 0.95) {
      print(paste0("% de données conservées : ", round(mean(keep), 4) * 100), collapse = "")
      print(paste0("% de probs conservées : ", round(sum(joint_density$joint_distr$p[keep]), 4) * 100), collapse = "")
    }
    
    if(any(!keep) & verbose) {
      print(paste0("nrow joint distr : ", nrow(joint_density$joint_distr), collapse = ""))
    } else if(sum(joint_density$joint_distr$p[keep]) < 0.95) {
      print(paste0("nrow joint distr : ", nrow(joint_density$joint_distr), collapse = ""))
    }
  } else {
    keep <- joint_density$joint_distr$p > 0
  }
  
  joint_density$joint_distr <- joint_density$joint_distr[keep, ]
  joint_density$grid_id <- joint_density$grid_id[keep, ]
  joint_density$grid_id <- as.data.frame(
    lapply(joint_density$grid_id, function(col) col - min(col) + 1)
  )
  joint_density$joint_distr$p <- joint_density$joint_distr$p / sum(joint_density$joint_distr$p)
  
  joint_density$domains <- lapply(joint_density$joint_distr[1:(ncol(joint_density$joint_distr) - 1)], function(x) sort(unique(x)))
  
  joint_density
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


drift_exact <- function(joint_density, clusters, a = 0.03) {
  priori <- new_priori(clusters)
  
  # arrondir pour accélérer les calculs en faisant du ==
  tmp <- joint_density$joint_distr
  tmp <- tmp[, -ncol(tmp), drop = FALSE]
  tmp <- round(tmp, 2)
  priori[, "mu"] <- round(priori[, "mu"], 2)
  
  priori_likely <- apply(tmp, 1, function(x) prod(priori[match(x, priori[, "mu"]), "p"]))
  
  res <- (1 - a) * joint_density$joint_distr$p + a * priori_likely
  joint_density$joint_distr$p <- res / sum(res)
  joint_density
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
  grid <- as.data.frame(new_priori(clusters))
  
  grid_id <- data.frame(1:nrow(grid))
  colnames(grid_id) <- name
  joint_density_init <- new_priori(clusters)[, "p"]
  joint_distr <- data.frame(grid["mu"], p = joint_density_init)
  colnames(joint_distr) <- c(name, "p")
  
  dom <- list(new_priori(clusters)[, "mu"])
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
          if (any(abs(distr[, "mu"] - x_i) < 0.01)) {
            distr[abs(distr[, "mu"] - x_i) < 0.01, "p"]
          } else {
            0
          }
        }))
      })),
    ncol = 2,
    dimnames = dimnames(init_distr())
  )
}

join_clusters <- function(clusters, which_i) {
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
  
  grid <- do.call(
    rbind,
    lapply(
      1:nrow(new_grid_id),
      function(i) {
        a <- new_grid_id[i, ]
        unlist(mapply(function(clust, j, id_col) clust$joint_distr[j, -id_col], clusters_subset, a, id_col_p, SIMPLIFY = FALSE))
      }
    )
  )
  
  new_grid_id_wide <- do.call(
    rbind,
    lapply(
      1:nrow(new_grid_id),
      function(i) {
        a <- new_grid_id[i, ]
        unlist(mapply(function(clust, j, id_col) clust$grid[j, ], clusters_subset, a, id_col_p, SIMPLIFY = FALSE))
      }
    )
  )
  
  prob_joint <- sapply(
    1:nrow(new_grid_id),
    function(i) {
      a <- new_grid_id[i, ]
      prod(mapply(function(clust, j, id_col) clust$joint_distr[j, id_col], clusters_subset, a, id_col_p))
    }
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

joint_density <- list(
  grid_id = data.frame(
    Louis = c(1, 1, 2),
    Alex = c(1, 2, 2),
    Jean = c(1, 1, 1)
  ),
  joint_distr = data.frame(
    Louis = c(10, 10, 20),
    Alex = c(10, 20, 20),
    Jean = c(10, 10, 10),
    p = c(0.3, 0.45, 0.25)
  ),
  domains = list(
    Louis = c(10, 20),
    Alex = c(10, 20),
    Jean = 10
  ),
  names = c("Louis", "Alex", "Jean")
)

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



recluster_dependancy <- function(joint_density, MI_thresh_for_indep = 0.05,
                                 max_cluster_size = Inf, min_cluster_size = 1) {
  if (length(joint_density$domains) == 1) return(list(joint_density))
  # faire un gros cluster et pour chq élément regarder si on peut les retirer d'une manière greedy
  combins <- generate_partitions(length(joint_density$domains), 2, min_cluster_size = min_cluster_size)
  MI <- sapply(combins, function(idx_out) 
    mutual_information(joint_density, idx_out)
  )
  if (length(MI) > 0) {
    
    
    # print(paste0("MI : ", paste0(round(MI, 3), collapse = ", ")))
    if (min(MI) < (MI_thresh_for_indep * log(length(joint_density$domains)))) {
      # l'exclure du gros cluster in
      # TODO moyen d'optimiser car dédouble les calculs déjà faits dans mutual_info()
      new_clusters <- marginal_joint_dependancy(joint_density, combins[[which.min(MI)]], format = 2)
      return(
        c(
          new_clusters[2],
          recluster_dependancy(new_clusters[[1]], MI_thresh_for_indep = MI_thresh_for_indep)
        )
      )
    }
    if (length(joint_density$domains) > max_cluster_size) {
      MI <- MI[length(joint_density$domains) - sapply(combins, length) <= max_cluster_size]
      combins <- combins[length(joint_density$domains) - sapply(combins, length) <= max_cluster_size]
      new_clusters <- marginal_joint_dependancy(
        joint_density,
        combins[[which.min(MI)]],
        format = 2
      )
      return(
        c(
          new_clusters[2],
          recluster_dependancy(
            new_clusters[[1]],
            MI_thresh_for_indep = MI_thresh_for_indep
          )
        )
      )
    }
  }
  joint_density$names <- unlist(names(joint_density$domains))
  list(joint_density)
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
  
  marginal_prob_in <- sapply(a_in, function(a) sum(joint_density$joint_distr$p[a == b_in]))
  marginal_prob_out <- sapply(a_out, function(a) sum(joint_density$joint_distr$p[a == b_out]))
  
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

mutual_information <- function(joint_density, idx_out) {
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
  
  deno_in <- marginal_prob_in[sapply(a_in, function(a) unlist(which(a == b_in)))]
  deno_out <- marginal_prob_out[sapply(a_out, function(a) unlist(which(a == b_out)))]
  
  sum(joint_density$joint_distr$p * (log(joint_density$joint_distr$p) - log(deno_in) - log(deno_out)))
}


compute_credibility <- function(distr, k = 0.05) {
  e <- sum(distr[, "mu"] * distr[, "p"])
  v <- (sum(distr[, "mu"]^2 * distr[, "p"]) - e^2)
  v <- max(v, 0) # car imprécisions sur les opérations floating point
  2 * pnorm(k * 50 / sqrt(v)) - 1 #2 * pnorm(k * e / sqrt(v)) - 1
}

is_exact_score_used_for_player <- function(distr, seuil = 0.7) {
  compute_credibility(distr) < seuil
}

