/*

  Extract genome sequences from a series of FASTA-formatted sequences.
  The sequence coordinates are taken from an input BED file.
  Optionally, the program computes and ouputs the nucleotide composition,
  in which case sequences can be extracted directly from the FASTA input or,
  as before, specified in a BED file.
  
  # Arguments:
  #   BED File [optional if -c option is set]
  #   Species  (e.g. hg19) [optional if -c optiom is set]
  # Options:
  #   -c Compute base composition [forward strand]
  #   -b Compute base composition for both strands  [-c mode set] 
  #   -r Compute base composition on reverse strand [-c mode set]

  Giovanna Ambrosini, EPFL/SV, giovanna.ambrosini@epfl.ch

  Copyright (c) 2015
  School of Life Sciences
  Ecole Polytechnique Federale de Lausanne
  and Swiss Institute of Bioinformatics
  EPFL SV ISREC UPNAE
  Station 15
  CH-1015 Lausanne, Switzerland.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

*/
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>
#include <ctype.h>
#include <limits.h>
#include "hashtable.h"
#ifdef DEBUG
#include <mcheck.h>
#endif

#define BUF_SIZE 8192
//#define BUF_SIZE 4194304 /* 4MB */
#define THIRTY_TWO_MEG 0x2000000ULL
#define FOUR_MEG 4000000
#define LINE_SIZE 1024
#define AC_MAX 18
#define CHR_NB 10
#define CHR_MAX 30
#define POS_MAX 16
#define INDEL_MAX 105
#define NUCL  5
#define LMAX  100
#define HDR_MAX 132
#define BED_RECORDS 1000
#define NB_OF_CHRS 40

typedef struct _options_t {
  int help;
  int debug;
  int nohdr;
  int bcomp;
  int both;
  int rev;
  int acPipe;
  char *dbPath;
} options_t;

static options_t options;

static char nucleotide[] = {'A','C','G','T', 'N'};
static char r_nucleotide[] = {'T','G','C','A', 'N'};

typedef struct _seq_t {
  char *hdr;
  int *seq;
  unsigned long len;
  char *ac;
} seq_t, *seq_p_t;

typedef struct _bed_t {
  unsigned long start;
  unsigned long end;
  char strand;
} bed_t, *bed_p_t;

typedef struct _chr_bed_t {
  bed_p_t bed_array;
} chr_bed_t, *chr_bed_p_t;

static chr_bed_t chr_record[NB_OF_CHRS];
static int bed_rec_cnt[NB_OF_CHRS] = {0};

FILE *fasta_in;
char *bedFile = NULL;

char *Species = NULL;

int TotChrom = 0;

static hash_table_t *ac_table = NULL;

static int
process_ac()
{
  FILE *input;
  int c;
  char buf[LINE_SIZE];
  char *chrFile;
  char chrom[12];
  int cLen;

  if (options.dbPath != NULL) {
      cLen = (int)strlen(options.dbPath) + strlen(Species) + 12;
      if ((chrFile = (char*)malloc(cLen * sizeof(char))) == NULL) {
        perror("process_ac: malloc");
        exit(1);
      }
      strcpy(chrFile, options.dbPath);
  } else {
      cLen = 21 + strlen(Species) + 12;
      if ((chrFile = (char*)malloc(cLen * sizeof(char))) == NULL) {
        perror("process_ac: malloc");
        exit(1);
      }
      strcpy(chrFile, "/home/local/db/genome");
  }
  strcat(chrFile, "/");
  strcat(chrFile, Species);
  strcat(chrFile, "/chr_NC_gi");

  input = fopen(chrFile, "r");
  if (input == NULL) {
    fprintf(stderr, "Could not open file %s: %s(%d)\n",
            chrFile, strerror(errno), errno);
    return 1;
  }

  do {
      c = fgetc(input);
  } while(c != '\n');

  ac_table = hash_table_new(MODE_COPY);
  while (fgets(buf, LINE_SIZE, input) != NULL) {
    char *s;
    char chr_nb[CHR_NB] = "";
    char ncbi_ac[AC_MAX] = "";
    int i = 0;
    int nb_len = 0;
    int ac_len = 0;
    s = buf;
    /* Check line */
    /* Get first character: if # skip line */
    if (*s == '#')
      continue;
    /* Chrom NB */
    while (*s != 0 && !isspace(*s)) {
      if (i >= CHR_NB) {
        fprintf(stderr, "AC too long in %s\n", s);
        fclose(input);
        exit(1);
      }
      chr_nb[i++] = *s++;
    }
    if (i < CHR_NB)
      chr_nb[i] = 0;
    nb_len = i + 1;
    while (isspace(*s))
      s++;
    /* Chromosome NCBI AC */
    i = 0;
    while (*s != 0 && !isspace(*s)) {
      if (i >= AC_MAX) {
        fprintf(stderr, "AC too long \"%s\" \n", s);
        fclose(input);
        exit(1);
      }
      ncbi_ac[i++] = *s++;
    }
    if (i < AC_MAX)
      ncbi_ac[i] = 0;
    ac_len = i + 1;
    strcpy(chrom, chr_nb);
    nb_len = (int)strlen(chrom) + 1;
    /* printf("adding key %s (len = %d) value %s (ac) (len = %d) to hash table\n", chrom, nb_len, ncbi_ac, ac_len); */
    /* Store both NCBI identifier to chrom number and chrom number to chrom number keys */
    hash_table_add(ac_table, ncbi_ac, (size_t)ac_len, chrom, (size_t)nb_len);
    if (options.debug) {
      char *cn = hash_table_lookup(ac_table, ncbi_ac, (size_t)ac_len);
      fprintf (stderr, " AC Hash table: %s (len = %d) -> %s (len = %d)\n", ncbi_ac, ac_len, cn, nb_len);
    }
    hash_table_add(ac_table, chrom, (size_t)nb_len, chrom, (size_t)nb_len);
    if (options.debug) {
      char *cn = hash_table_lookup(ac_table, chrom, (size_t)nb_len);
      fprintf (stderr, " AC Hash table: %s (len = %d) -> %s (len = %d)\n", chrom, nb_len, cn, nb_len);
    }
  }
  return 0;
}

static void 
change_chrnb(char *chrnb)
{
  if ( (strcmp(Species, "hg18") == 0) || (strcmp(Species, "hg19") == 0) || (strcmp(Species, "hg38") == 0) ) {
    if (strcmp(chrnb, "X") == 0)
      strcpy(chrnb, "23"); 
    if (strcmp(chrnb, "Y") == 0)
      strcpy(chrnb, "24"); 
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "25"); 
  } else if ( (strcmp(Species, "mm8") == 0) || (strcmp(Species, "mm9") == 0) || (strcmp(Species, "mm10") == 0) ) {
    if (strcmp(chrnb, "X") == 0)
      strcpy(chrnb, "20"); 
    if (strcmp(chrnb, "Y") == 0)
      strcpy(chrnb, "21"); 
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "22"); 
  } else if ( (strcmp(Species, "bosTau3") == 0) || (strcmp(Species, "bosTau8") == 0) ) {
    if (strcmp(chrnb, "X") == 0)
      strcpy(chrnb, "30"); 
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "31"); 
  } else if ( (strcmp(Species, "canFam2") == 0) || (strcmp(Species, "canFam3") == 0) ) {
    if (strcmp(chrnb, "X") == 0)
      strcpy(chrnb, "39"); 
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "40"); 
  } else if ( (strcmp(Species, "panTro2") == 0) || (strcmp(Species, "panTro5") == 0) ) {
    if (strcmp(chrnb, "2A") == 0)
      strcpy(chrnb, "2"); 
    if (strcmp(chrnb, "2B") == 0)
      strcpy(chrnb, "3"); 
    if (strcmp(chrnb, "3") == 0)
      strcpy(chrnb, "23"); 
    if (strcmp(chrnb, "X") == 0)
      strcpy(chrnb, "24"); 
    if (strcmp(chrnb, "Y") == 0)
      strcpy(chrnb, "25"); 
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "26"); 
  } else if ( (strcmp(Species, "rn5") == 0) ) {
    if (strcmp(chrnb, "X") == 0)
      strcpy(chrnb, "21"); 
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "22"); 
  } else if ( (strcmp(Species, "rn6") == 0) ) {
    if (strcmp(chrnb, "X") == 0)
      strcpy(chrnb, "21"); 
    if (strcmp(chrnb, "Y") == 0)
      strcpy(chrnb, "22"); 
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "23"); 
  } else if ( (strcmp(Species, "amel5") == 0) ) {
    if (strcmp(chrnb, "LG1") == 0)
      strcpy(chrnb, "1"); 
    if (strcmp(chrnb, "LG2") == 0)
      strcpy(chrnb, "2"); 
    if (strcmp(chrnb, "LG3") == 0)
      strcpy(chrnb, "3"); 
    if (strcmp(chrnb, "LG4") == 0)
      strcpy(chrnb, "4"); 
    if (strcmp(chrnb, "LG5") == 0)
      strcpy(chrnb, "5"); 
    if (strcmp(chrnb, "LG6") == 0)
      strcpy(chrnb, "6"); 
    if (strcmp(chrnb, "LG7") == 0)
      strcpy(chrnb, "7"); 
    if (strcmp(chrnb, "LG8-24") == 0)
      strcpy(chrnb, "8"); 
    if (strcmp(chrnb, "LG9") == 0)
      strcpy(chrnb, "9"); 
    if (strcmp(chrnb, "LG10") == 0)
      strcpy(chrnb, "10"); 
    if (strcmp(chrnb, "LG11") == 0)
      strcpy(chrnb, "11"); 
    if (strcmp(chrnb, "LG12") == 0)
      strcpy(chrnb, "12"); 
    if (strcmp(chrnb, "LG13") == 0)
      strcpy(chrnb, "13"); 
    if (strcmp(chrnb, "LG14") == 0)
      strcpy(chrnb, "14"); 
    if (strcmp(chrnb, "LG15") == 0)
      strcpy(chrnb, "15"); 
    if (strcmp(chrnb, "LG16") == 0)
      strcpy(chrnb, "16"); 
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "17"); 
  } else if ( (strcmp(Species, "dm3") == 0) ) {
    if (strcmp(chrnb, "2L") == 0)
      strcpy(chrnb, "1"); 
    if (strcmp(chrnb, "2R") == 0)
      strcpy(chrnb, "2"); 
    if (strcmp(chrnb, "3L") == 0)
      strcpy(chrnb, "3"); 
    if (strcmp(chrnb, "3R") == 0)
      strcpy(chrnb, "4"); 
    if (strcmp(chrnb, "4") == 0)
      strcpy(chrnb, "5"); 
    if (strcmp(chrnb, "X") == 0)
      strcpy(chrnb, "6"); 
  } else if ( (strcmp(Species, "dm6") == 0) ) {
    if (strcmp(chrnb, "2L") == 0)
      strcpy(chrnb, "1"); 
    if (strcmp(chrnb, "2R") == 0)
      strcpy(chrnb, "2"); 
    if (strcmp(chrnb, "3L") == 0)
      strcpy(chrnb, "3"); 
    if (strcmp(chrnb, "3R") == 0)
      strcpy(chrnb, "4"); 
    if (strcmp(chrnb, "4") == 0)
      strcpy(chrnb, "5"); 
    if (strcmp(chrnb, "X") == 0)
      strcpy(chrnb, "6"); 
    if (strcmp(chrnb, "Y") == 0)
      strcpy(chrnb, "7"); 
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "8"); 
  } else if ( (strcmp(Species, "danRer7") == 0) || (strcmp(Species, "danRer10") == 0) ) {
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "26"); 
  } else if ( (strcmp(Species, "susScr3") == 0) ) {
    if (strcmp(chrnb, "X") == 0)
      strcpy(chrnb, "19"); 
    if (strcmp(chrnb, "Y") == 0)
      strcpy(chrnb, "20"); 
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "21"); 
  } else if ( (strcmp(Species, "ce6") == 0) || (strcmp(Species, "ce10") == 0) || (strcmp(Species, "ce11") == 0) ) {
    if (strcmp(chrnb, "I") == 0)
      strcpy(chrnb, "1"); 
    if (strcmp(chrnb, "II") == 0)
      strcpy(chrnb, "2"); 
    if (strcmp(chrnb, "III") == 0)
      strcpy(chrnb, "3"); 
    if (strcmp(chrnb, "IV") == 0)
      strcpy(chrnb, "4"); 
    if (strcmp(chrnb, "V") == 0)
      strcpy(chrnb, "5"); 
    if (strcmp(chrnb, "X") == 0)
      strcpy(chrnb, "6"); 
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "7"); 
  } else if ( (strcmp(Species, "spo2") == 0) ) {
    if (strcmp(chrnb, "I") == 0)
      strcpy(chrnb, "1"); 
    if (strcmp(chrnb, "II") == 0)
      strcpy(chrnb, "2"); 
    if (strcmp(chrnb, "III") == 0)
      strcpy(chrnb, "3"); 
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "4"); 
  } else if ( (strcmp(Species, "oryLat") == 0) ) {
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "25"); 
  } else if ( (strcmp(Species, "oreNil2") == 0) ) {
    if (strcmp(chrnb, "LG1") == 0)
      strcpy(chrnb, "1"); 
    if (strcmp(chrnb, "LG2") == 0)
      strcpy(chrnb, "2"); 
    if (strcmp(chrnb, "LG3") == 0)
      strcpy(chrnb, "3"); 
    if (strcmp(chrnb, "LG4") == 0)
      strcpy(chrnb, "4"); 
    if (strcmp(chrnb, "LG5") == 0)
      strcpy(chrnb, "5"); 
    if (strcmp(chrnb, "LG6") == 0)
      strcpy(chrnb, "6"); 
    if (strcmp(chrnb, "LG7") == 0)
      strcpy(chrnb, "7"); 
    if (strcmp(chrnb, "LG8-24") == 0)
      strcpy(chrnb, "8"); 
    if (strcmp(chrnb, "LG9") == 0)
      strcpy(chrnb, "9"); 
    if (strcmp(chrnb, "LG10") == 0)
      strcpy(chrnb, "10"); 
    if (strcmp(chrnb, "LG11") == 0)
      strcpy(chrnb, "11"); 
    if (strcmp(chrnb, "LG12") == 0)
      strcpy(chrnb, "12"); 
    if (strcmp(chrnb, "LG13") == 0)
      strcpy(chrnb, "13"); 
    if (strcmp(chrnb, "LG14") == 0)
      strcpy(chrnb, "14"); 
    if (strcmp(chrnb, "LG15") == 0)
      strcpy(chrnb, "15"); 
    if (strcmp(chrnb, "LG16-21") == 0)
      strcpy(chrnb, "16"); 
    if (strcmp(chrnb, "LG17") == 0)
      strcpy(chrnb, "17"); 
    if (strcmp(chrnb, "LG18") == 0)
      strcpy(chrnb, "18"); 
    if (strcmp(chrnb, "LG19") == 0)
      strcpy(chrnb, "19"); 
    if (strcmp(chrnb, "LG20") == 0)
      strcpy(chrnb, "20"); 
    if (strcmp(chrnb, "LG21") == 0)
      strcpy(chrnb, "21"); 
    if (strcmp(chrnb, "LG22") == 0)
      strcpy(chrnb, "22"); 
    if (strcmp(chrnb, "LG23") == 0)
      strcpy(chrnb, "23"); 
    if (strcmp(chrnb, "MT") == 0)
      strcpy(chrnb, "24"); 
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "24"); 
  } else if ( (strcmp(Species, "xenTro9") == 0) ) {
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "11"); 
  } else if ( (strcmp(Species, "sacCer2") == 0) || (strcmp(Species, "sacCer3") == 0) ) {
    if (strcmp(chrnb, "I") == 0)
      strcpy(chrnb, "1"); 
    if (strcmp(chrnb, "II") == 0)
      strcpy(chrnb, "2"); 
    if (strcmp(chrnb, "III") == 0)
      strcpy(chrnb, "3"); 
    if (strcmp(chrnb, "IV") == 0)
      strcpy(chrnb, "4"); 
    if (strcmp(chrnb, "V") == 0)
      strcpy(chrnb, "5"); 
    if (strcmp(chrnb, "VI") == 0)
      strcpy(chrnb, "6"); 
    if (strcmp(chrnb, "VII") == 0)
      strcpy(chrnb, "7"); 
    if (strcmp(chrnb, "VIII") == 0)
      strcpy(chrnb, "8"); 
    if (strcmp(chrnb, "IX") == 0)
      strcpy(chrnb, "9"); 
    if (strcmp(chrnb, "X") == 0)
      strcpy(chrnb, "10"); 
    if (strcmp(chrnb, "XI") == 0)
      strcpy(chrnb, "11"); 
    if (strcmp(chrnb, "XII") == 0)
      strcpy(chrnb, "12"); 
    if (strcmp(chrnb, "XIII") == 0)
      strcpy(chrnb, "13"); 
    if (strcmp(chrnb, "XIV") == 0)
      strcpy(chrnb, "14"); 
    if (strcmp(chrnb, "XV") == 0)
      strcpy(chrnb, "15"); 
    if (strcmp(chrnb, "XVI") == 0)
      strcpy(chrnb, "16"); 
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "17"); 
  } else if ( (strcmp(Species, "araTha1") == 0) ) {
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "6"); 
  } else if ( (strcmp(Species, "orySat") == 0) ) {
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "13"); 
    if (strcmp(chrnb, "Pltd") == 0)
      strcpy(chrnb, "14"); 
    if (strcmp(chrnb, "B1") == 0)
      strcpy(chrnb, "15"); 
  } else if ( (strcmp(Species, "zm3") == 0) ) {
    if (strcmp(chrnb, "M") == 0)
      strcpy(chrnb, "11"); 
    if (strcmp(chrnb, "Pltd") == 0)
      strcpy(chrnb, "12"); 
  } else {
    fprintf (stderr, "Invalid Species\n");
    exit(1);
  }
}

static void
load_bed(const char *file)
{
  char buf[LINE_SIZE];
  FILE *f;
  int i;
  size_t mLen[NB_OF_CHRS];

  if (options.debug)
    fprintf (stderr, "Prosessing BED File %s...\n", file);

  if ((f = fopen(file, "r")) == NULL) {
    perror("load_bed: fopen");
    exit(1);
  }
  /* Allocate SNP Records for each chromosome */
  for (i = 0; i < NB_OF_CHRS; i++) {
    mLen[i] = BED_RECORDS;
    if ( (chr_record[i].bed_array = (bed_p_t) malloc(mLen[i] * sizeof(bed_t))) == NULL ) {
      perror("load_bed: malloc chr_record.bed_array");
      exit(1);
    }
  }
  while (fgets(buf, LINE_SIZE, f) != NULL) {
    char *s;
    char chr_nb[CHR_NB] = "";
    char start_p[POS_MAX] = "";
    char end_p[POS_MAX] = "";
    char strand = '\0';
    s = buf;
    i = 0;
    /* Chrom NB */
    //printf ("BED Line: %s\n", s);
    while (*s != 0 && !isspace(*s)) {
      if (i >= CHR_NB) {
        fprintf(stderr, "CHR NB too long in %s\n", s);
        fclose(f);
        exit(1);
      }
      if (i > 2)
        chr_nb[i-3] = *s;
      i++;
      s++;
    }
    if (i < CHR_NB)
      chr_nb[i-3] = 0;
    //printf ("%s: Chr NB: %s\n", Species, chr_nb);
    change_chrnb(chr_nb);
    while (isspace(*s))
      s++;
    /* Start Position */
    i = 0;
    while (*s != 0 && !isspace(*s)) {
      if (i >= POS_MAX) {
        fprintf(stderr, "RS too long in %s\n", s);
        fclose(f);
        exit(1);
      }
      start_p[i++] = *s++;
    }
    if (i < POS_MAX)
      start_p[i] = 0;
    //printf ("START POS: %s\n", rs_id);
    while (isspace(*s))
      s++;
    /* End Position */
    i = 0;
    while (*s != 0 && !isspace(*s)) {
      if (i >= POS_MAX) {
        fprintf(stderr, "Position too long in %s\n", s);
        fclose(f);
        exit(1);
      }
      end_p[i++] = *s++;
    }
    if (i < POS_MAX)
      end_p[i] = 0;
    //printf ("END POS: %s\n", pos);
    while (isspace(*s))
      s++;
    /* Skip Name Field  */
    while (*s != 0 && !isspace(*s))
      s++;
    while (isspace(*s))
      s++;
    /* Skip Score Field */
    while (*s != 0 && !isspace(*s))
      s++;
    while (isspace(*s))
      s++;
    /* Strand           */
    strand = *s;

    int chr = atoi(chr_nb);
    if (chr-1 >= 0) {
      if ((size_t)bed_rec_cnt[chr-1] >= mLen[chr-1]) {
        mLen[chr-1] *= 2;
        if ( (chr_record[chr-1].bed_array = (bed_p_t) realloc(chr_record[chr-1].bed_array, mLen[chr-1] * sizeof(bed_t))) == NULL ) {
          perror("load_bed: realloc chr_record.bed_array");
          exit(1);
        }
      }
      //printf ("Copying START POS %d END POS %d STRAND %c to BED REC #%d of CHROM NB %d\n", atoi(start_p), atoi(end_p), strand, bed_rec_cnt[chr-1], chr); 
      chr_record[chr-1].bed_array[bed_rec_cnt[chr-1]].start = atoi(start_p);
      chr_record[chr-1].bed_array[bed_rec_cnt[chr-1]].end = atoi(end_p);
      chr_record[chr-1].bed_array[bed_rec_cnt[chr-1]].strand = strand;
      //printf ("Chr %d RECORD (#%d): %lu  %lu  %c\n", chr, chr_record[chr-1].bed_array[bed_rec_cnt[chr-1]].start, chr_record[chr-1].bed_array[bed_rec_cnt[chr-1]].end, strand);
      bed_rec_cnt[chr-1]++;
    }
  }
  fclose(f);
}

void 
dump_bed()
{
  for (int i = 1; i <= NB_OF_CHRS; i++) {
    //printf ("Chromosome %d: # of BED Records %d\n", i, bed_rec_cnt[i-1]);
    for (int k = 0; k < bed_rec_cnt[i-1]; k++) {
      fprintf (stderr, "%d\t%lu\t%lu\t%c\n", i, chr_record[i-1].bed_array[k].start, chr_record[i-1].bed_array[k].end, chr_record[i-1].bed_array[k].strand);
    }
  }
}

static int
process_seqs(FILE *input, const char *iFile)
{
  char buf[BUF_SIZE], *res;
  seq_t seq;
  size_t mLen;
  char *chr_nb;
  int ac_len;
  int chr;

  if (input == NULL) {
    FILE *f = fopen(iFile, "r");
    if (f == NULL) {
      fprintf(stderr, "Could not open file %s: %s(%d)\n",
  	    iFile, strerror(errno), errno);
      return -1;
    }
    input = f;
  }
  if (options.debug != 0)
    fprintf(stderr, "Processing file %s\n", iFile);

  while ((res = fgets(buf, BUF_SIZE, input)) != NULL
	 && buf[0] != '>')
    ;
  if (res == NULL || buf[0] != '>') {
    fprintf(stderr, "Could not find a sequence in file %s\n", iFile);
    if (input != stdin) {
      fclose(input);
    }
    return 1;
  }
  seq.hdr = malloc(HDR_MAX * sizeof(char));
  seq.ac = malloc(AC_MAX * sizeof(char));
  seq.seq = malloc(THIRTY_TWO_MEG * sizeof(int));
  mLen = THIRTY_TWO_MEG;
  while (res != NULL) {
    /* Get the header */
    if (buf[0] != '>') {
      fprintf(stderr, "Could not find a sequence header in file %s\n", iFile);
      if (input != stdin) {
        fclose(input);
      }
      return 1;
    }
    char *s = buf;
    s += 1;
    int i = 0;
    while (*s && !isspace(*s)) {
      if (i >= HDR_MAX) {
        fprintf(stderr, "Fasta Header too long \"%s\" in file %s\n", res, iFile);
        fclose(input);
        return 1;
      }
      seq.hdr[i++] = *s++;
    }
    if (i < HDR_MAX)
      seq.hdr[i] = 0;
    /* Get AC  */
    s = seq.hdr;
    for (i = 0; i < options.acPipe; i++) {
      s = strchr(s, '|');
      if (s == NULL) {
        fprintf(stderr, "Bad header line \"%s\" in file %s\n", res, iFile);
        fclose(input);
        return 1;
      }
      s += 1;
    }
    if (options.acPipe == 0)
      s = seq.hdr;
    seq.len = 0;
    while (*s && *s != '|' && *s != ';' && !isspace(*s)) {
      if (seq.len >= AC_MAX) {
        fprintf(stderr, "AC too long \"%s\" in file %s\n", res, iFile);
        fclose(input);
        return 1;
      }
      seq.ac[seq.len++] = *s++;
    }
    if (seq.len < AC_MAX)
      seq.ac[seq.len] = 0;
    /* Gobble sequence  */ 
    seq.len = 0;
    while ((res = fgets(buf, BUF_SIZE, input)) != NULL && buf[0] != '>') {
      char c;
      int n;
      s = buf;
      while ((c = *s++) != 0) {
        if (isalpha(c)) {
          c = (char) toupper(c);
          switch (c) {
          case 'A':
            n = 0;
            break;
          case 'C':
            n = 1;
            break;
          case 'G':
            n = 2;
            break;
          case 'T':
            n = 3;
            break;
          case 'N':
            n = 4;
            break;
          default:
            n = 4;
            ;
          }
          if (seq.len >= mLen) {
            mLen *= 2;
            seq.seq = realloc(seq.seq, (size_t)mLen * sizeof(int));
            if (seq.seq == NULL) {
              perror("process_seqs: realloc");
              exit(1);
            }
          }
          seq.seq[seq.len++] = n;
        }
      }
    }
    /* We now have the (not nul terminated) sequence.
       Process it. */
    if (seq.len != 0) {
      /* Process BED file  */
      /* Get Chromosome number */
      ac_len = (int)strlen(seq.ac) + 1; 
      chr_nb = hash_table_lookup(ac_table, seq.ac, (size_t)ac_len);
      if (chr_nb == NULL)
        continue;
      change_chrnb(chr_nb);
      chr = atoi(chr_nb);
      if (options.debug != 0) {
        fprintf (stderr, "Processing BED file for Chromosome %d\n", chr);
        fprintf (stderr, "Number of BED Records %d\n", bed_rec_cnt[chr-1]);
      }
      /* Loop on BED Record Array for Chromosome chr   */
      for (int k = 0; k < bed_rec_cnt[chr-1]; k++) {
        unsigned long start = chr_record[chr-1].bed_array[k].start;
        unsigned long end = chr_record[chr-1].bed_array[k].end;
        // Print Sequence Header
        printf(">%s [%lu..%lu]\n", seq.hdr, start, end);
        int cnt = 0;
        if (chr_record[chr-1].bed_array[k].strand == '-') {
          for (unsigned int i = end-1; i >= start-1; i--) {
            cnt++;
            printf("%c", r_nucleotide[seq.seq[i]]);
            if ( ((cnt) % 70) == 0 ) {
              printf("\n");
            }
          }
        } else {
          for (unsigned int i = start-1; i < end; i++) {
            cnt++;
            printf("%c", nucleotide[seq.seq[i]]);
            if ( ((cnt) % 70) == 0 ) {
              printf("\n");
            }
          }
        }
        printf("\n");
      } /* End loop on BED Records  */
    }   /* If Seq Length not NULL   */
  }
  free(seq.seq);
  fclose(input);
  return 0;
}

static int
compute_bcomp_r(FILE *input, const char *iFile)
{
  char buf[BUF_SIZE], *res;
  seq_t seq;
  size_t mLen;
  char *chr_nb;
  int ac_len;
  int chr;
  unsigned int bcomp[5] = {0, 0, 0, 0, 0};
  unsigned long tot_len = 0;

  if (input == NULL) {
    FILE *f = fopen(iFile, "r");
    if (f == NULL) {
      fprintf(stderr, "Could not open file %s: %s(%d)\n",
  	    iFile, strerror(errno), errno);
      return -1;
    }
    input = f;
  }
  if (options.debug != 0)
    fprintf(stderr, "Processing file %s\n", iFile);

  while ((res = fgets(buf, BUF_SIZE, input)) != NULL
	 && buf[0] != '>')
    ;
  if (res == NULL || buf[0] != '>') {
    fprintf(stderr, "Could not find a sequence in file %s\n", iFile);
    if (input != stdin) {
      fclose(input);
    }
    return 1;
  }
  seq.hdr = malloc(HDR_MAX * sizeof(char));
  seq.ac = malloc(AC_MAX * sizeof(char));
  seq.seq = malloc(THIRTY_TWO_MEG * sizeof(int));
  mLen = THIRTY_TWO_MEG;
  while (res != NULL) {
    /* Get the header */
    if (buf[0] != '>') {
      fprintf(stderr, "Could not find a sequence header in file %s\n", iFile);
      if (input != stdin) {
        fclose(input);
      }
      return 1;
    }
    char *s = buf;
    s += 1;
    int i = 0;
    while (*s && !isspace(*s)) {
      if (i >= HDR_MAX) {
        fprintf(stderr, "Fasta Header too long \"%s\" in file %s\n", res, iFile);
        fclose(input);
        return 1;
      }
      seq.hdr[i++] = *s++;
    }
    if (i < HDR_MAX)
      seq.hdr[i] = 0;
    /* Get AC  */
    s = seq.hdr;
    for (i = 0; i < options.acPipe; i++) {
      s = strchr(s, '|');
      if (s == NULL) {
        fprintf(stderr, "Bad header line \"%s\" in file %s\n", res, iFile);
        fclose(input);
        return 1;
      }
      s += 1;
    }
    if (options.acPipe == 0)
      s = seq.hdr;
    seq.len = 0;
    while (*s && *s != '|' && *s != ';' && !isspace(*s)) {
      if (seq.len >= AC_MAX) {
        fprintf(stderr, "AC too long \"%s\" in file %s\n", res, iFile);
        fclose(input);
        return 1;
      }
      seq.ac[seq.len++] = *s++;
    }
    if (seq.len < AC_MAX)
      seq.ac[seq.len] = 0;
    /* Gobble sequence  */ 
    seq.len = 0;
    while ((res = fgets(buf, BUF_SIZE, input)) != NULL && buf[0] != '>') {
      char c;
      int n;
      s = buf;
      while ((c = *s++) != 0) {
        if (isalpha(c)) {
          c = (char) toupper(c);
          switch (c) {
          case 'A':
            n = 0;
            break;
          case 'C':
            n = 1;
            break;
          case 'G':
            n = 2;
            break;
          case 'T':
            n = 3;
            break;
          case 'N':
            n = 4;
            break;
          default:
            n = 4;
            ;
          }
          if (seq.len >= mLen) {
            mLen *= 2;
            seq.seq = realloc(seq.seq, (size_t)mLen * sizeof(int));
            if (seq.seq == NULL) {
              perror("process_seqs: realloc");
              exit(1);
            }
          }
          seq.seq[seq.len++] = n;
        }
      }
    }
    /* We now have the (not nul terminated) sequence.
       Process it. */
    if (seq.len != 0) {
      if (bedFile != NULL) {
        /* Process BED file  */
        /* Get Chromosome number */
        ac_len = (int)strlen(seq.ac) + 1; 
        chr_nb = hash_table_lookup(ac_table, seq.ac, (size_t)ac_len);
        if (chr_nb == NULL)
          continue;
        change_chrnb(chr_nb);
        chr = atoi(chr_nb);
        if (options.debug != 0) {
          fprintf (stderr, "Processing BED file for Chromosome %d\n", chr);
          fprintf (stderr, "Number of BED Records %d\n", bed_rec_cnt[chr-1]);
        }
        /* Loop on BED Record Array for Chromosome chr   */
        for (int k = 0; k < bed_rec_cnt[chr-1]; k++) {
          unsigned long start = chr_record[chr-1].bed_array[k].start;
          unsigned long end = chr_record[chr-1].bed_array[k].end;
          // Print Sequence Header
          if (chr_record[chr-1].bed_array[k].strand == '-') {
            for (unsigned int i = start-1; i < end; i++) {
              bcomp[seq.seq[i]]++;
              tot_len++;
            }
          } else {
            for (unsigned int i = start -1; i < end; i++) {
              if (seq.seq[i] < 4)
                bcomp[3-seq.seq[i]]++;
              else
                bcomp[seq.seq[i]]++;
              tot_len++;
            }
          }
        } /* End loop on BED Records  */
      } else {  /* Process the entire sequence  */
        //printf("Seq lenght: %lu\n", seq.len);
        for (unsigned int i = 0; i < seq.len; i++) {
          if (seq.seq[i] < 4)
            bcomp[3-seq.seq[i]]++; 
          else
            bcomp[seq.seq[i]]++; 
        }
        tot_len +=seq.len;
      }
    }   /* If Seq length not NULL   */
  }
  //printf("A:%d , C:%d , G:%d , T:%d , N:%d\n", bcomp[0], bcomp[1], bcomp[2], bcomp[3], bcomp[4]);
  fprintf(stderr, "Total Sequence length: %lu\n", tot_len);
  printf("%.4f,%.4f,%.4f,%.4f\n", (double)((double)(bcomp[0]+bcomp[4]/4)/tot_len), (double)((double)(bcomp[1]+bcomp[4]/4)/tot_len), (double)((double)(bcomp[2]+bcomp[4]/4)/tot_len), (double)((double)(bcomp[3]+bcomp[4]/4)/tot_len)); 
  free(seq.seq);
  fclose(input);
  return 0;
}

static int
compute_bcomp(FILE *input, const char *iFile)
{
  char buf[BUF_SIZE], *res;
  seq_t seq;
  size_t mLen;
  char *chr_nb;
  int ac_len;
  int chr;
  unsigned int bcomp[5] = {0, 0, 0, 0, 0};
  unsigned long tot_len = 0;

  if (input == NULL) {
    FILE *f = fopen(iFile, "r");
    if (f == NULL) {
      fprintf(stderr, "Could not open file %s: %s(%d)\n",
  	    iFile, strerror(errno), errno);
      return -1;
    }
    input = f;
  }
  if (options.debug != 0)
    fprintf(stderr, "Processing file %s\n", iFile);

  while ((res = fgets(buf, BUF_SIZE, input)) != NULL
	 && buf[0] != '>')
    ;
  if (res == NULL || buf[0] != '>') {
    fprintf(stderr, "Could not find a sequence in file %s\n", iFile);
    if (input != stdin) {
      fclose(input);
    }
    return 1;
  }
  seq.hdr = malloc(HDR_MAX * sizeof(char));
  seq.ac = malloc(AC_MAX * sizeof(char));
  seq.seq = malloc(THIRTY_TWO_MEG * sizeof(int));
  mLen = THIRTY_TWO_MEG;
  while (res != NULL) {
    /* Get the header */
    if (buf[0] != '>') {
      fprintf(stderr, "Could not find a sequence header in file %s\n", iFile);
      if (input != stdin) {
        fclose(input);
      }
      return 1;
    }
    char *s = buf;
    s += 1;
    int i = 0;
    while (*s && !isspace(*s)) {
      if (i >= HDR_MAX) {
        fprintf(stderr, "Fasta Header too long \"%s\" in file %s\n", res, iFile);
        fclose(input);
        return 1;
      }
      seq.hdr[i++] = *s++;
    }
    if (i < HDR_MAX)
      seq.hdr[i] = 0;
    /* Get AC  */
    s = seq.hdr;
    for (i = 0; i < options.acPipe; i++) {
      s = strchr(s, '|');
      if (s == NULL) {
        fprintf(stderr, "Bad header line \"%s\" in file %s\n", res, iFile);
        fclose(input);
        return 1;
      }
      s += 1;
    }
    if (options.acPipe == 0)
      s = seq.hdr;
    seq.len = 0;
    while (*s && *s != '|' && *s != ';' && !isspace(*s)) {
      if (seq.len >= AC_MAX) {
        fprintf(stderr, "AC too long \"%s\" in file %s\n", res, iFile);
        fclose(input);
        return 1;
      }
      seq.ac[seq.len++] = *s++;
    }
    if (seq.len < AC_MAX)
      seq.ac[seq.len] = 0;
    /* Gobble sequence  */ 
    seq.len = 0;
    while ((res = fgets(buf, BUF_SIZE, input)) != NULL && buf[0] != '>') {
      char c;
      int n;
      s = buf;
      while ((c = *s++) != 0) {
        if (isalpha(c)) {
          c = (char) toupper(c);
          switch (c) {
          case 'A':
            n = 0;
            break;
          case 'C':
            n = 1;
            break;
          case 'G':
            n = 2;
            break;
          case 'T':
            n = 3;
            break;
          case 'N':
            n = 4;
            break;
          default:
            n = 4;
            ;
          }
          if (seq.len >= mLen) {
            mLen *= 2;
            seq.seq = realloc(seq.seq, (size_t)mLen * sizeof(int));
            if (seq.seq == NULL) {
              perror("process_seqs: realloc");
              exit(1);
            }
          }
          seq.seq[seq.len++] = n;
        }
      }
    }
    /* We now have the (not nul terminated) sequence.
       Process it. */
    if (seq.len != 0) {
      if (bedFile != NULL) {
        /* Process BED file  */
        /* Get Chromosome number */
        ac_len = (int)strlen(seq.ac) + 1; 
        chr_nb = hash_table_lookup(ac_table, seq.ac, (size_t)ac_len);
        if (chr_nb == NULL)
          continue;
        change_chrnb(chr_nb);
        chr = atoi(chr_nb);
        if (options.debug != 0) {
          fprintf (stderr, "Processing BED file for Chromosome %d\n", chr);
          fprintf (stderr, "Number of BED Records %d\n", bed_rec_cnt[chr-1]);
        }
        /* Loop on BED Record Array for Chromosome chr   */
        for (int k = 0; k < bed_rec_cnt[chr-1]; k++) {
          unsigned long start = chr_record[chr-1].bed_array[k].start;
          unsigned long end = chr_record[chr-1].bed_array[k].end;
          // Print Sequence Header
          if (chr_record[chr-1].bed_array[k].strand == '-') {
            for (unsigned int i = start-1; i < end; i++) {
              if (seq.seq[i] < 4)
                bcomp[3-seq.seq[i]]++; 
              else
                bcomp[seq.seq[i]]++;
              tot_len++;
            }
          } else {
            for (unsigned int i = start-1; i < end; i++) {
              bcomp[seq.seq[i]]++;
              tot_len++;
            }
          }
        } /* End loop on BED Records  */
      } else {  /* Process the entire sequence  */
        //printf("Seq lenght: %lu\n", seq.len);
        for (unsigned int i = 0; i < seq.len; i++) {
          bcomp[seq.seq[i]]++;
        }
        tot_len +=seq.len;
      }
    }   /* If Seq length not NULL   */
  }
  //printf("A:%d , C:%d , G:%d , T:%d , N:%d\n", bcomp[0], bcomp[1], bcomp[2], bcomp[3], bcomp[4]);
  fprintf(stderr, "Total Sequence length: %lu\n", tot_len);
  if (options.both) {
    double bcomp_at = (double)((double)(bcomp[0]+bcomp[4]/4)/tot_len);
    double bcomp_cg = (double) 0.5 - bcomp_at;
    printf("%.2f,%.2f,%.2f,%.2f\n", bcomp_at, bcomp_cg, bcomp_cg, bcomp_at); 
  } else {
    printf("%.4f,%.4f,%.4f,%.4f\n", (double)((double)(bcomp[0]+bcomp[4]/4)/tot_len), (double)((double)(bcomp[1]+bcomp[4]/4)/tot_len), (double)((double)(bcomp[2]+bcomp[4]/4)/tot_len), (double)((double)(bcomp[3]+bcomp[4]/4)/tot_len)); 
  }
  free(seq.seq);
  fclose(input);
  return 0;
}

int
main(int argc, char *argv[])
{
#ifdef DEBUG
  mcheck(NULL);
  mtrace();
#endif
  options.acPipe = 2;
  options.dbPath = NULL;
  while (1) {
    int c = getopt(argc, argv, "dhbcri:f:p:s:");
    if (c == -1)
      break;
    switch (c) {
    case 'd':
      options.debug = 1;
      break;
    case 'h':
      options.help = 1;
      break;
    case 'b':
      options.both = 1;
      break;
    case 'c':
      options.bcomp = 1;
      break;
    case 'i':
      options.acPipe = atoi(optarg);
      break;
    case 'f':
      bedFile = optarg;
      break;
    case 'p':
      options.dbPath = optarg;
      break;
    case 'r':
      options.rev = 1;
      break;
    case 's':
      Species = optarg;
      break;
    case '?':
      break;
    default:
      printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }
  if (optind > argc || options.help == 1 || (bedFile == NULL && options.bcomp == 0) || (Species == NULL && options.bcomp == 0)) {
    fprintf(stderr,
	    "Usage: %s [options] [-f <bed_file>] [-s <species>] [<] [<fasta_file>|stdin]\n"
	    "      where options are:\n"
	    "        -d          Print debug information\n"
	    "        -h          Show this help text\n"
            "        -c          Compute base composition [def=on forward strand]\n"
            "        -b          Compute base composition on both strand [-c is required]\n"
            "        -r          Compute base composition on reverse strand [-c is required]\n"
            "        -i <int>    AC index (after how many pipes |) for FASTA header [%d]\n"
            "        -p <path>   Use <path> to locate the chr_NC_gi file [if BED file is given]\n"
            "                    [default is: $HOME/db/genome]\n"
	    "\n\tExtract BED regions from a set of FASTA-formatted sequences.\n"
            "\tThe extracted sequences are written to standard output.\n"
            "\tOptionally (-c), the program computes and only outputs the base composition,\n"
            "\tin which case sequences can be extracted directly from the FASTA input or,\n"
            "\tas for the extraction mode, specified in a BED file. If the BED file is not\n"
            "\tgiven, the <species> argument is not required. If base composition mode is set\n"
            "\t(-c option), the program can optionally compute it on both strands (-b option)\n"
            "\tfor strand-symmetric base composition or on the reverse strand only (-r).\n\n",
	    argv[0], options.acPipe);
    return 1;
  }
  if (argc > optind) {
      if(!strcmp(argv[optind],"-")) {
          fasta_in = stdin;
      } else {
          fasta_in = fopen(argv[optind], "r");
          if (fasta_in == NULL) {
              fprintf(stderr, "Unable to open '%s': %s(%d)\n",
                  argv[optind], strerror(errno), errno);
             exit(EXIT_FAILURE);
          }
          if (options.debug)
             fprintf(stderr, "Processing file %s\n", argv[optind]);
      }
  } else {
      fasta_in = stdin;
  }
  if (options.debug != 0) {
    if (fasta_in != stdin) {
      fprintf(stderr, "FASTA sequence file : %s\n", argv[optind]);
    } else {
      fprintf(stderr, "FASTA sequence file from STDIN\n");
    }
    if (bedFile != NULL)
      fprintf(stderr, "Extract sequences from BED file %s\n", bedFile);
    else
      fprintf(stderr, "Process the entire FASTA sequence file\n");
    fprintf(stderr, "Species Assembly %s\n", Species);
  }
  if (bedFile != NULL) {
    load_bed(bedFile);
    if (options.debug)
      dump_bed();
    if (process_ac() == 0) {
      if (options.debug)
        fprintf(stderr, " HASH Table for chromosome access identifier initialized\n");
    } else {
      return 1;
    }
  }
  if (options.bcomp) { 
    if (options.rev) {
      if (compute_bcomp_r(fasta_in, argv[optind++]) != 0)
        return 1;
    } else {
      if (compute_bcomp(fasta_in, argv[optind++]) != 0)
        return 1;
    }
  } else {
    if (process_seqs(fasta_in, argv[optind++]) != 0)
      return 1;
  }
  return 0;
}
