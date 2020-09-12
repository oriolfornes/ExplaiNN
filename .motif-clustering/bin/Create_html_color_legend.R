################
## Read files ##
################
# taxa <- c("insects", "vertebrates", "nematodes", "fungi", "plants")
taxa <- c("vertebrates")

sapply(taxa, function(taxon){
  
  ## Read the table with the motif information (Id, link, color, TF class)
  table.files.dir <- file.path("/home/jamondra/Downloads/Tmp/jaspar_2020/")
  Jaspar.table.file <- paste("JASPAR2016_", taxon, "_matrices_for_matrix_clustering.tsv", sep = "")
  Jaspar.table.file <- "jaspar2020_nematodes_annotations_with_colours.tsv"
  Jaspar.table.file <- file.path(table.files.dir, Jaspar.table.file)
  
  Jaspar.table <-read.csv(Jaspar.table.file, sep = "\t", header = FALSE)
  colnames(Jaspar.table) <- c("ID", "Link", "Color", "Class")
  
  ## Create an association table Color-TF class
  Color.Class.table <- unique(Jaspar.table[,c(3,4)])
  
  # Color.Class.table <- Color.Class.table[order(Color.Class.table$Color),]
  
  #######################
  ## Export HTML table ##
  #######################
  
  ## Table header + tail
  head.tab <- "<div id='Color_class_tab' style='display: inline-block;float:left;position:relative;' class='color-legend' width='375px'><p style='font-size:12px;padding:0px;border:0px'><b></b></p><table id='Color_class_table' class='hover compact stripe' cellspacing='0' width='375px' style='padding:15px;align:center;'><thead><tr><th > Color </th><th> TF Class </th> </tr></thead><tbody>"
  tab.lines <- paste("\n<tr><td class='color-box' style='background-color: --color--';></td><td>--TFClass--</td></tr>", collapse = "")
  tail.tab <- "</tbody></table></div>"
  
  ## Table body
  table.body <- sapply(1:nrow(Color.Class.table), function(r.nb){
    
    ## Set variables
    tab.lines.current.line <- tab.lines
    TF.class.Color <- Color.Class.table[r.nb,1]
    TF.class.Name <- Color.Class.table[r.nb,2]
    
    tab.lines.current.line <- gsub("--TFClass--", TF.class.Name, tab.lines.current.line)
    tab.lines.current.line <- gsub("--color--", TF.class.Color, tab.lines.current.line)
    
  })
  table.body <- paste(table.body, collapse = "")
  
  html.table <- paste(head.tab, table.body, tail.tab, collapse = "")
  writeLines(html.table)
  
})