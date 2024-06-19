source("src/import.R")
source("src/update_scores.R")
source("src/plots.R")


#show_current_ranking(update_scores(players, scores[1:34, ]))

#TODO games played per players
#TODO points gained/loss per players

tmp <- show_ranking_history(scores)
players <- tmp[[1]]
ggplot(tmp[[2]])+
  geom_line(aes(x=date, y=score, col=player), linewidth=1)+
  geom_point(aes(x=date, y=score, col=player))+
  scale_color_discrete(breaks = tmp[[2]][order(tmp[[2]][tmp[[2]][, 1] == max(tmp[[2]][, 1]), 3], decreasing = T), 2])+
  theme_bw()

show_current_ranking(players)

show_current_probs(players)

show_skill_level(players)

show_detailed_skill(players)

lapply(show_skill_level(players) / 100, transition_matrix)

