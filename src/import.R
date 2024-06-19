source("src/players.R")
source("src/scores.R")

scores <- scores_init()

players <- list()

data <- read.csv2("data/data.csv")
colnames(data) <- c("joueur_A1", "joueur_A2", "joueur_B1", "joueur_B2", "win", "score_A", "score_B", "game_len", "date")
data <- data[, 1:9]
data <- data[-1, ]
data <- data[data[, "joueur_B1"] != "", ]

mapping_joueurs <- matrix(c("W", "Will",
                            "É", "Éti",
                            "V", "Vic",
                            "P", "Phil",
                            "A", "Ant",
                            "X", "Xav",
                            "G", "Gab",
                            "CS", "Clau",
                            "J", "Jon",
                            "CC", "Charl",
                            "N", "Nate",
                            "Z", "Zach",
                            "LG", "Lor",
                            "JT", "Jas"), ncol=2, byrow=T)

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
  
  scores <- add_scores(
    subdata[, c(1:5, 8)],
    scores,
    date = d
  )
}

