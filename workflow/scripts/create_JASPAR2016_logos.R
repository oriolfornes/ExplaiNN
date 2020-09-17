library(TFBSTools);
library(JASPAR2016);

opts <- list();
#opts[["species"]] <- 9606;
#opts[["name"]] <- "RUNX1";
#opts[["type"]] <- "SELEX";
#opts[["all_versions"]] <- TRUE;
opts[["all"]] <- TRUE;
PFMatrixList = getMatrixSet(JASPAR2016, opts);

library(ggseqlogo);

svg_dir <-  "./JASPAR2014_logos_SVG/";
png_small_dir <- "./JASPAR2014_logos_PNG_small/";
png_large_dir <- "./JASPAR2014_logos_PNG_large/";
dir.create(svg_dir);
dir.create(png_small_dir);
dir.create(png_large_dir);
for(indx in sequence(length(PFMatrixList))){
    mat <- PFMatrixList[indx];
    name <- ID(mat);
    svg(paste0(svg_dir, name, ".svg"), width=18.66);
    p <- ggseqlogo(Matrix(mat), font="roboto_regular");
    print(p);
    null <- dev.off();
    svg(paste0(svg_dir, name, ".rc.svg"), width=18.66);
    p <- ggseqlogo(reverseComplement(Matrix(mat)[[name]]), font="roboto_regular");
    print(p);
    null <- dev.off();
    png(paste0(png_small_dir, name, ".png"), height=100);
    p <- ggseqlogo(Matrix(mat)[[name]], "font=roboto_regular");
    print(p);
    null <- dev.off();
    png(paste0(png_small_dir, name, ".rc.png"), height=100);
    p <- ggseqlogo(reverseComplement(Matrix(mat)[[name]]), font="roboto_regular");
    print(p);
    null <- dev.off();

    png(paste0(png_large_dir, name, ".png"), height=180, res=150);
    p <- ggseqlogo(Matrix(mat)[[name]], "font=roboto_regular");
    print(p);
    null <- dev.off();
    png(paste0(png_large_dir, name, ".rc.png"), height=180, res=150);
    p <- ggseqlogo(reverseComplement(Matrix(mat)[[name]]), font="roboto_regular");
    print(p);
    null <- dev.off();
}
