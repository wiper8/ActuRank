include_exact_points <- F

source("src/import.R")
source("src/update_scores.R")
source("src/plots.R")
source("src/recommend.R")

scores_stats(scores, players)

games_matchups(scores, players)

tmp <- show_ranking_history(scores)
players <- tmp[[1]]

#ranking history
ggplot(tmp[[2]])+
  geom_line(aes(x=date, y=score, col=player), linewidth=1)+
  geom_point(aes(x=date, y=score, col=player))+
  scale_color_discrete(breaks = tmp[[2]][order(tmp[[2]][tmp[[2]][, 1] == max(tmp[[2]][, 1]), 3], decreasing = T), 2])+
  theme_bw()

sort(sapply(players, compute_credibility), decreasing = T)

show_current_ranking(players)

show_detailed_skill(players)

show_skill_level(players)

recommend_next_game(players, names_present = NULL)

show_current_probs(players)

lapply(show_skill_level(players) / 100, transition_matrix)

show_played_against_grid(players, scores)

