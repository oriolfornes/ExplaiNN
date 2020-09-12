#####################
## Load R packages ##
#####################
required.libraries <- c("data.table",
                        "dplyr")

for (lib in required.libraries) {
  if (!require(lib, character.only=TRUE)) {
    install.packages(lib)
    suppressPackageStartupMessages(library(lib, character.only=TRUE))
  }
}


## Read motif information matrix
motif.info.tab <- read.csv("/home/jamondra/Downloads/All_motifs_dm3.tab", sep = "\t", header = FALSE)
colnames(motif.info.tab) <- c("Experiment", "TF", "TF2", "Centrimo_file", "Centrality_pval", "Logo")

motif.info.tab$Centrality_pval <- as.numeric(as.vector(motif.info.tab$Centrality_pval))
motif.info.tab.summary <- motif.info.tab %>% 
                            filter(Centrality_pval <= -200) %>% 
                            group_by(TF) %>% 
                            tally(sort = TRUE)


View(
motif.info.tab %>% 
  filter(Centrality_pval <= -100) %>% 
  group_by(TF) %>% 
  tally(sort = TRUE)
)