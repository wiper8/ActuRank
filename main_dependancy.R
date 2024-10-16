dataset <- "ping"

if (dataset == "ping") {
  include_exact_points <- FALSE
  dim_len_mu <- 18
}
if (dataset == "spike") {
  include_exact_points <- TRUE
  dim_len_mu <- 17
}
if (dataset == "pickle") {
  include_exact_points <- TRUE
  dim_len_mu <- 18
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

games_matchups(scores, players)
set.seed(2024L)

# generate_GIF_images(scores)

tmp <- show_ranking_history_dependancy(scores, dataset)
players <- tmp[[1]]
clusters <- tmp[[3]]
#TODO retirer Frédéric

# score history
ggplot(tmp[[2]])+
  theme_bw()+
  geom_line(aes(x=date, y=score, col=player), linewidth=1)+
  geom_point(aes(x=date, y=score, col=player), data=tmp[[2]][tmp[[2]]$played, ])+
  scale_color_discrete(breaks = tmp[[2]][order(tmp[[2]][tmp[[2]][, 1] == max(tmp[[2]][, 1]), "score"], decreasing = T), "player"])+
  coord_cartesian(xlim = c(min(as.Date(scores$date)), max(as.Date(scores$date))))

ggplot(tmp[[2]])+
  theme_bw()+
  geom_line(aes(x=date, y=skill, col=player), linewidth=1)+
  geom_point(aes(x=date, y=skill, col=player), data=tmp[[2]][tmp[[2]]$played, ])+
  scale_color_discrete(breaks = tmp[[2]][order(tmp[[2]][tmp[[2]][, 1] == max(tmp[[2]][, 1]), "score"], decreasing = T), "player"])+
  coord_cartesian(xlim = c(min(as.Date(scores$date)), max(as.Date(scores$date))))

# ranking history
ggplot(tmp[[2]])+
  theme_bw()+
  geom_line(aes(x=day_i, y=rank, col=player), linewidth=1)+
  geom_point(aes(x=day_i, y=rank, col=player), data=tmp[[2]][tmp[[2]]$played, ])+
  scale_color_discrete(breaks = tmp[[2]][order(tmp[[2]][tmp[[2]][, 1] == max(tmp[[2]][, 1]), "score"], decreasing = T), "player"])+
  theme(panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = 1:max(tmp[[2]][, "day_i"], na.rm = TRUE))+
  scale_y_reverse(breaks = 1:max(tmp[[2]][, "rank"], na.rm = TRUE))

ggplot(tmp[[2]])+
  theme_bw()+
  geom_line(aes(x=day_i, y=rank_skill, col=player), linewidth=1)+
  geom_point(aes(x=day_i, y=rank_skill, col=player), data=tmp[[2]][tmp[[2]]$played, ])+
  scale_color_discrete(breaks = tmp[[2]][order(tmp[[2]][tmp[[2]][, 1] == max(tmp[[2]][, 1]), "score"], decreasing = T), "player"])+
  theme(panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = 1:max(tmp[[2]][, "day_i"], na.rm = TRUE))+
  scale_y_reverse(breaks = 1:max(tmp[[2]][, "rank_skill"], na.rm = TRUE))

sort(sapply(players, compute_credibility), decreasing = T)

show_current_ranking(clusters, scores)

show_current_ranking(clusters, scores, show_credibility = TRUE)

show_detailed_skill(players)

show_detailed_skill_per_player(players)

show_IC_skill(players)

show_skill_level(players)

recommend_next_game(players, names_present = NULL)

show_current_probs_exact(clusters)

lapply(show_skill_level(players), transition_matrix)

show_played_against_grid(players, scores)

