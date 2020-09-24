BEGIN{
    OFS = "\t";
}

!/^#/{
    split($1, chr, ":");
    split(chr[2], pos, "-");
    if ($4 == "D"){
        strand = "+";
    }else{
        strand = "-";
    }
    start = pos[1] + $5 - 1;
    end = pos[1] + $6;
    name = genome"_"chr[1]":"start+1"-"end"("strand")";
    print chr[1], start, end, name, ".", strand;
}