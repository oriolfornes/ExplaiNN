BEGIN{
    OFS = "\t";
    if(! species){
        species="hg38";
    }
}

!/^#/{
    split($1, chr, ":");
    split(chr[2], pos, "-");
    if ($4 == "D"){
        strand = "+";
    }else{
        if ($4 == "R"){
            strand = "-";
        }
    }
    start = pos[1] + $5 - 1;
    end = pos[1] + $6;
    print chr[1], start, end, species"_"chr[1]":"start+1"-"end"("strand")", ".", strand;
}
