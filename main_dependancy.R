dataset <- "spike"
granularity_level <- 1

ping_granularity <- data.frame(
  level = c(1, 2, 3, 4),
  # TODO
  dim_len_mu = c(10, 20, 30, 40),
  # TODO
  seuil_freq = c(0, 20, 50, 300)
)
if (dataset == "spike") {
  include_exact_points <- TRUE
  dim_len_mu <- 20
}
if (dataset == "pickle") {
  include_exact_points <- TRUE
  dim_len_mu <- 30
}

if (dataset == "ping") {
  source("src/import_ping.R")
  # retirer vieux data
  scores <- scores[scores$date >= as.Date("2024-01-01"), ]
  #TODO essayer dajouter 10 wins de xav contre vic pour voir cque ca fait
} else if (dataset == "spike") {
  source("src/import_spike.R")
} else if (dataset == "pickle") {
  source("src/import_pickle.R")
  source("src/pickle.R")
}
if (!"serve_for_pt" %in% colnames(scores)) {
  scores$serve_for_pt <- FALSE
}
source("src/update_scores.R")
source("src/plots.R")
source("src/recommend.R")

scores_stats(scores, players)
nrow(scores)
sum(scores$score_A + scores$score_B)

games_matchups(scores, players)
set.seed(2024L)

if (dataset == "ping") {
  include_exact_points <- FALSE
  dim_len_mu <- ping_granularity[granularity_level, "dim_len_mu"]
  seuil_freq <- ping_granularity[granularity_level, "seuil_freq"]
  
  names_peu_freq <- scores_stats(scores, players)$games_played
  names_peu_freq <- names(names_peu_freq[names_peu_freq <= seuil_freq])
  scores <- scores[apply(scores[, 2:5], 1, function(x) !any(x %in% names_peu_freq)), ]
}

# generate_GIF_images(scores)
tmp <- show_ranking_history_dependancy(scores, dataset)
players <- tmp[[1]]
graph_data <- tmp[[2]]
clusters <- tmp[[3]]

# score history
ggplot(graph_data)+
  theme_bw()+
  geom_line(aes(x=date, y=score, col=player), linewidth=1)+
  geom_point(aes(x=date, y=score, col=player), data=graph_data[graph_data$played, ])+
  scale_color_discrete(breaks = graph_data[order(graph_data[graph_data[, 1] == max(graph_data[, 1]), "score"], decreasing = T), "player"])+
  coord_cartesian(xlim = c(min(as.Date(scores$date)), max(as.Date(scores$date))))

ggplot(graph_data)+
  theme_bw()+
  geom_line(aes(x=date, y=skill, col=player), linewidth=1)+
  geom_point(aes(x=date, y=skill, col=player), data=graph_data[graph_data$played, ])+
  scale_color_discrete(breaks = graph_data[order(graph_data[graph_data[, 1] == max(graph_data[, 1]), "score"], decreasing = T), "player"])+
  coord_cartesian(xlim = c(min(as.Date(scores$date)), max(as.Date(scores$date))))

# ranking history
ggplot(graph_data)+
  theme_bw()+
  geom_line(aes(x=day_i, y=rank, col=player), linewidth=1)+
  geom_point(aes(x=day_i, y=rank, col=player), data=graph_data[graph_data$played, ])+
  scale_color_discrete(breaks = graph_data[order(graph_data[graph_data[, 1] == max(graph_data[, 1]), "score"], decreasing = T), "player"])+
  theme(panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = 1:max(graph_data[, "day_i"], na.rm = TRUE))+
  scale_y_reverse(breaks = 1:max(graph_data[, "rank"], na.rm = TRUE))

ggplot(graph_data)+
  theme_bw()+
  geom_line(aes(x=day_i, y=rank_skill, col=player), linewidth=1)+
  geom_point(aes(x=day_i, y=rank_skill, col=player), data=graph_data[graph_data$played, ])+
  scale_color_discrete(breaks = graph_data[order(graph_data[graph_data[, 1] == max(graph_data[, 1]), "score"], decreasing = T), "player"])+
  theme(panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = 1:max(graph_data[, "day_i"], na.rm = TRUE))+
  scale_y_reverse(breaks = 1:max(graph_data[, "rank_skill"], na.rm = TRUE))

ggplot(graph_data)+
  theme_bw()+
  geom_line(aes(x=day_i, y=credibility, col=player), linewidth=1)+
  geom_point(aes(x=day_i, y=credibility, col=player), data=graph_data[graph_data$played, ])+
  scale_color_discrete(breaks = graph_data[order(graph_data[graph_data[, 1] == max(graph_data[, 1]), "score"], decreasing = T), "player"])+
  theme(panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = 1:max(graph_data[, "day_i"], na.rm = TRUE))+
  ylim(0, 1)

sort(
  sapply(
    clusters,
    function(distr) compute_multivariate_credibility(distr$joint_distr)),
  decreasing = TRUE
)

show_current_ranking(clusters, scores)

show_current_ranking(clusters, scores, show_credibility = TRUE)

show_detailed_skill(players)

show_detailed_skill_per_player(players)

test_hyp(clusters)

show_IC_skill(players)

show_skill_level(players)

recommend_next_game(players, names_present = NULL)

show_current_probs_exact2(clusters)

lapply(show_skill_level(players), transition_matrix)

show_played_against_grid(players, scores)

