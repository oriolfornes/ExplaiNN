#####################
## Load R packages ##
#####################
required.libraries <- c("data.table",
                        "dplyr",
                        "optparse",
                        "RColorBrewer")

for (lib in required.libraries) {
  if (!require(lib, character.only=TRUE)) {
    install.packages(lib)
    suppressPackageStartupMessages(library(lib, character.only=TRUE))
  }
}


####################
## Read arguments ##
####################
option_list = list(
  make_option(c("-i", "--in_jaspar_table"), type="character", default=NULL, 
              help="Input JASPAR table", metavar="character"),
  make_option(c("-c", "--colour_html"), type="character", default=NULL, 
              help="Colour HTML table", metavar="character"),
  make_option(c("-o", "--out_jaspar_table_colours"), type="character", default=NULL, 
              help="Output JASPAR table with TF class-colour assignation", metavar="character")
  
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

jaspar.tab.file <- opt$in_jaspar_table
jaspar.tab.colour.file <- opt$out_jaspar_table_colours
jaspar.html.tab.colour <- opt$colour_html


## Debug
# jaspar.tab.file <- "/home/jamondra/Downloads/Tmp/jaspar_2020/jaspar2020_fungi_annotations.tsv"
# jaspar.tab.file <- "/home/jamondra/Downloads/Tmp/jaspar_2020/jaspar2020_insects_annotations.tsv"
# jaspar.tab.file <- "/home/jamondra/Downloads/Tmp/jaspar_2020/jaspar2020_nematodes_annotations.tsv"
# jaspar.tab.file <- "/home/jamondra/Downloads/Tmp/jaspar_2020/jaspar2020_plants_annotations.tsv"
# jaspar.tab.file <- "/home/jamondra/Downloads/Tmp/jaspar_2020/jaspar2020_vertebrates_annotations.tsv"
# jaspar.tab.colour.file <- "/home/jamondra/Downloads/Tmp/jaspar_2020/jaspar2020_vertebrates_annotations_with_colors.tsv"

#######################
## Read JASPAR table ##
#######################
message("; Reading JASPAR table")
jaspar.tab <- fread(jaspar.tab.file)
jaspar.tab <- within(jaspar.tab, rm("colour"))

## Replace "." by "_" in the matrix ID
jaspar.tab$matrix_id <- gsub(jaspar.tab$matrix_id, pattern = "\\.", replacement = "_")

## For the dimers considers only the first TF class
jaspar.tab$class <- gsub(jaspar.tab$class, pattern = ",.+$", replacement = "", perl = T)
jaspar.tab$class <- gsub(jaspar.tab$class, pattern = "::.+$", replacement = "", perl = T)

# table(jaspar.tab$class)
# View(jaspar.tab)

##################################
## TF class - Color assignation ##
##################################
message("; Assigning colours to TF classes")
TF.classes <- as.vector(unique(jaspar.tab$class))
TF.classes.nb <- length(TF.classes)

nb.seed.colors <- ifelse(TF.classes.nb < 8,
                         yes = TF.classes.nb,
                         no = 8)
  
## Dark2 palette seed max = 8
class.colors <- colorRampPalette(brewer.pal(nb.seed.colors, "Dark2"), space="Lab")(TF.classes.nb)
names(class.colors) <- TF.classes

df.class.colour <- data.frame(colour = class.colors, class = names(class.colors))


##########################################
## Create HTML (TF class::Colour) table ##
##########################################

## Table header + tail
head.tab <- "<div id='Color_class_tab' style='display: inline-block;float:left;position:relative;' class='color-legend' width='450px'><p style='font-size:12px;padding:0px;border:0px'><b></b></p><table id='Color_class_table' class='hover compact stripe' cellspacing='0' width='450px' style='padding:15px;align:center;'><thead><tr><th > Color </th><th> TF Class </th> </tr></thead><tbody>"
tab.lines <- paste("\n<tr><td class='color-box' style='background-color: --color--';></td><td>--TFClass--</td></tr>", collapse = "")
tail.tab <- "<tr><td class='non_validated'>*</td><td>Non-validated</td></tr></tbody></table></div>"

## Table body
table.body <- sapply(1:nrow(df.class.colour), function(r.nb){
  
  ## Set variables
  tab.lines.current.line <- tab.lines
  TF.class.Color <- df.class.colour[r.nb,1]
  TF.class.Name <- df.class.colour[r.nb,2]
  
  tab.lines.current.line <- gsub("--TFClass--", TF.class.Name, tab.lines.current.line)
  tab.lines.current.line <- gsub("--color--", TF.class.Color, tab.lines.current.line)
  
})
table.body <- paste(table.body, collapse = "")

html.table <- paste(head.tab, table.body, tail.tab, collapse = "")
writeLines(html.table, jaspar.html.tab.colour)


##################
## Export table ##
##################
message("; Export output table: ", jaspar.tab.colour.file)

## Combine tables
jaspar.tab.colour <- merge(jaspar.tab, df.class.colour, by = "class")

## Order columns
jaspar.tab.colour <- jaspar.tab.colour[,c("matrix_id", "URL", "colour", "class", "name")]
fwrite(jaspar.tab.colour, file = jaspar.tab.colour.file, sep = "\t", quote = FALSE)

# cd /workspace/rsat/tmp/www-data/2019/08/27/matrix-clustering_2019-08-27.110825_UxjKx0
# sudo rsync -rptuvl rsat@rsat-tagc.univ-mrs.fr:/home/jcastro/Jaspar_2020/results/JASPAR_2020_matrix_clustering/* .