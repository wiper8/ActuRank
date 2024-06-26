source("src/players.R")
source("src/scores.R")


scores <- scores_init()

players <- list()

data <- read.csv2("data/data.csv")
colnames(data) <- c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2", "win", "score_A", "score_B", "game_len", "date")
data <- data[, 1:9]
data <- data[-1, ]
data <- data[data[, "joueur_B1"] != "", ]
data$score_A <- as.numeric(data$score_A)
data$score_B <- as.numeric(data$score_B)
data$win <- as.numeric(data$win)

mapping_joueurs <- matrix(c("W", "Will",
                            "É", "Éti",
                            "V", "Vic",
                            "P", "Phil",
                            "ANTP", "Ant",
                            "X", "Xav",
                            "G", "Gab",
                            "CS", "Clau",
                            "JOC", "Jon",
                            "CC", "Charl",
                            "N", "Nate",
                            "Z", "Zach",
                            "LG", "Lor",
                            "JT", "Jas",
                            "APA2", "AlexP",
                            "JCH", "Jacob",
                            "AUJ", "Audrey",
                            "CV", "Carl",
                            "M", "Mariève",
                            "LL", "Louis",
                            "JPL", "JPL",
                            "AR", "AlexRich"), ncol=2, byrow=T)

for(i in 1:4)
  data[, i] <- mapping_joueurs[match(data[, i], mapping_joueurs), 2]

noms <- unique(unlist(data[, 1:4]))
noms <- noms[noms != "" & !is.na(noms)]

for(n in noms) {
  players <- add_player(n, players)
}

#names(players) <- mapping_joueurs[match(c("W", "P", "É"), mapping_joueurs), 2]

data[, "date"] <- as.Date(data[, "date"], tryFormats = c("%d/%m/%Y", "%Y-%m-%d"))
dates <- sort(unique(data[, "date"]))

for(d in as.character(dates)) {
  print(d)
  subdata <- as.data.frame(data[data[, "date"] == d, ])
  
  #if(include_exact_points) {
  #  scores <- add_scores(
  #    do.call(rbind, apply(subdata[, 1:8], 1, scores_per_pt_converter)),
  #    scores,
  #    date = d
  #  )
  #} else {
    scores <- add_scores(
      subdata[, 1:8],
      scores,
      date = d
    )
  #}
  
}

name <- unique(unlist(scores[, c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2")]))
name <- name[!is.na(name)]

players_very_low_exposure <- sapply(name, function(n) {
  sum(as.numeric(scores[apply(scores, 1, function(score) {
    n %in% score[c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2")]
  }), "game_len"]))
})
players_very_low_exposure <- names(players_very_low_exposure)[players_very_low_exposure < 100]

