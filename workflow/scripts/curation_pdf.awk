BEGIN{
    FS = "\t";
    print "\\documentclass{article}";
    print "\\usepackage{graphicx}";
    print "\\usepackage{url}";
    print "\\pagestyle{empty}";
    print "\\begin{document}";
}

/^GSE/{
    repo = "20170705_GEO_MACS_Jeanne_localrsat/GEO_peakMax_perExp";
}

/^ENCSR/{
    repo = "20170704_ENCODE_MACS_Marius_localrsat/ENCODE_peakMax_perExp";
}

!/^PWM/{
    centrimo = $1".pssm.501bp.fa.sites.centrimo";
    print "\\section*{"$4"}";
    print "\\includegraphics{{"repo"/"$1"}.png}";
    print "\\includegraphics[scale=.2]{{"repo"/"centrimo"}.pdf}";
}

END{
    print "\\end{document}";
}
