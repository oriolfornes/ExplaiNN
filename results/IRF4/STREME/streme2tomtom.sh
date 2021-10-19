#!/usr/bin/env bash

tomtom -oc WT.tomtom WT/streme.txt JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
tomtom -oc T95R.tomtom T95R/streme.txt JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
tomtom -oc WT.201bp.tomtom WT.201bp/streme.txt JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
tomtom -oc T95R.201bp.tomtom T95R.201bp/streme.txt JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
tomtom -oc WT-unique.tomtom WT-unique/streme.txt JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
tomtom -oc T95R-unique.tomtom T95R-unique/streme.txt JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
tomtom -oc T95R-WT-intersect.tomtom T95R-WT-intersect/streme.txt JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
