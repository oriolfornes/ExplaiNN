#!/usr/bin/env bash

if ! [ -f hg38/hg38.fa.sizes ]; then
	genomepy install -p UCSC -g ./ -r "^chr[\dXYM]{1,2}$" -t 8 -f hg38
fi
