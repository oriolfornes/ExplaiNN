#!/usr/bin/env bash

FASTA_DIR=

tomtom -no-ssc -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 10.0 -oc WT.tomtom WT/streme.txt db/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme db/MOUSE/uniprobe_mouse.meme db/EUKARYOTE/jolma2013.meme
tomtom -no-ssc -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 10.0 -oc T95R.tomtom T95R/streme.txt db/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme db/MOUSE/uniprobe_mouse.meme db/EUKARYOTE/jolma2013.meme
tomtom -no-ssc -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 10.0 -oc WT.201bp.tomtom WT.201bp/streme.txt db/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme db/MOUSE/uniprobe_mouse.meme db/EUKARYOTE/jolma2013.meme
tomtom -no-ssc -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 10.0 -oc T95R.201bp.tomtom T95R.201bp/streme.txt db/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme db/MOUSE/uniprobe_mouse.meme db/EUKARYOTE/jolma2013.meme
tomtom -no-ssc -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 10.0 -oc WT-unique.tomtom WT-unique/streme.txt db/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme db/MOUSE/uniprobe_mouse.meme db/EUKARYOTE/jolma2013.meme
tomtom -no-ssc -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 10.0 -oc T95R-unique.tomtom T95R-unique/streme.txt db/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme db/MOUSE/uniprobe_mouse.meme db/EUKARYOTE/jolma2013.meme
tomtom -no-ssc -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 10.0 -oc T95R-WT-intersect.tomtom T95R-WT-intersect/streme.txt db/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme db/MOUSE/uniprobe_mouse.meme db/EUKARYOTE/jolma2013.meme
