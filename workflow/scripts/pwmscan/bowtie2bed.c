/*
  bowtie2bed.c

  Convert Bowtie output into BED format
  
  # Arguments:
  # score file
  # species
  # Bowtie output file 

  Giovanna Ambrosini, EPFL, Giovanna.Ambrosini@epfl.ch

  Copyright (c) 2014 EPFL and Swiss Institute of Bioinformatics.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
/*
#define DEBUG 
*/
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <ctype.h>
#include "hashtable.h"
#ifdef DEBUG
#include <mcheck.h>
#endif

/*#define BUF_SIZE 4096 */
#define BUF_SIZE 8192
#define LINE_SIZE 1024
#define CHR_NB 10
#define TAG_MAX 128
#define HDR_MAX 36
#define POS_MAX 16
#define EXT_MAX 128
#define SCORE_MAX 12

typedef struct _options_t {
  char *scFile;
  char *dbPath;
  char *colSep;
  int score;
  int norm;
  int mism;
  int db;
  int debug;
  int help;
} options_t;

static options_t options;

char *Species = NULL;
static hash_table_t *ac_table = NULL;
int tagLen = 0;
int misMatch = 0;

void 
reverse(char *str)
{
  char temp;
  int i = 0;
  int j = 0;
 
  j = strlen(str) - 1;
  while (i < j) {
    temp = str[i];
    str[i] = str[j];
    str[j] = temp;
    i++;
    j--;
  }
}

void 
complement(char *str)
{
  int i = 0;
  int j = 0;
 
  j = strlen(str);
  while (i < j) {
    switch (str[i]) {
      case 'A':
        str[i] = 'T';
        break;
      case 'B':
        str[i] = 'V';
        break;
      case 'C':
        str[i] = 'G';
        break;
      case 'D':
        str[i] = 'H';
        break;
      case 'G':
        str[i] = 'C';
        break;
      case 'H':
        str[i] = 'D';
        break;
      case 'M':
        str[i] = 'K';
        break;
      case 'N':
        str[i] = 'N';
        break;
      case 'R':
        str[i] = 'Y';
        break;
      case 'S':
        str[i] = 'S';
        break;
      case 'T':
        str[i] = 'A';
        break;
      case 'U':
        str[i] = 'A';
        break;
      case 'V':
        str[i] = 'B';
        break;
      case 'W':
        str[i] = 'W';
        break;
      case 'X':
        str[i] = 'X';
        break;
      case 'Y':
        str[i] = 'R';
        break;
      case 'a':
        str[i] = 't';
        break;
      case 'b':
        str[i] = 'v';
        break;
      case 'c':
        str[i] = 'g';
        break;
      case 'd':
        str[i] = 'h';
        break;
      case 'g':
        str[i] = 'c';
        break;
      case 'h':
        str[i] = 'd';
        break;
      case 'm':
        str[i] = 'k';
        break;
      case 'n':
        str[i] = 'n';
        break;
      case 'r':
        str[i] = 'y';
        break;
      case 's':
        str[i] = 's';
        break;
      case 't':
        str[i] = 'a';
        break;
      case 'u':
        str[i] = 'a';
        break;
      case 'v':
        str[i] = 'b';
        break;
      case 'w':
        str[i] = 'w';
        break;
      case 'x':
        str[i] = 'x';
        break;
      case 'y':
        str[i] = 'r';
        break;
      default:
        ;
    }
    i++;
  }
}

int 
process_ac() 
{
  FILE *input;
  int c;
  char buf[LINE_SIZE];
  char *chrFile;
  int cLen;

  if (options.db) {
      cLen = strlen(options.dbPath) + strlen(Species) + 12;
      if ((chrFile = malloc(cLen * sizeof(char))) == NULL) {
        perror("process_ac: malloc");
        exit(1);
      }
      strcpy(chrFile, options.dbPath);
  } else {
      cLen = 21 + strlen(Species) + 12;
      if ((chrFile = malloc(cLen * sizeof(char))) == NULL) {
        perror("process_ac: malloc");
        exit(1);
      }
      strcpy(chrFile, "/home/local/db/genome");
  }
  strcat(chrFile, "/"); 
  strcat(chrFile, Species); 
  strcat(chrFile, "/chr_hdr"); 

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
  //while (fscanf(input, "%s\t%s\t%d\n", chr_nb, ncbi_ac, &gi) != EOF) {
  while (fgets(buf, LINE_SIZE, input) != NULL) {
    char *s;
    char chr_nb[CHR_NB] = "";
    char ncbi_hdr[HDR_MAX] = "";
    int i = 0;
    int nb_len = 0;
    int hdr_len = 0;
    /* Chrom NB */
    s = buf;
    while (*s != 0 && !isspace(*s)) {
      if (i >= CHR_NB) {
        fprintf(stderr, "Chr NB too long in %s\n", s);
        fclose(input);
        exit(1);
      }
      chr_nb[i++] = *s++;
    }
    if (i < CHR_NB)
      chr_nb[i] = 0;
    //nb_len = strlen(chr_nb) + 1;
    nb_len = i + 1;
    while (isspace(*s))
    s++; 
    /* Chromosome AC */ 
    i = 0;
    while (*s != 0 && !isspace(*s)) {
      if (i >= HDR_MAX) {
        fprintf(stderr, "Seq Header too long \"%s\" \n", s);
        fclose(input);
        exit(1);
      }
      ncbi_hdr[i++] = *s++;
    }
    if (i < HDR_MAX)
      ncbi_hdr[i] = 0;
    //ac_len = strlen(ncbi_ac) + 1;
    hdr_len = i + 1;
    hash_table_add(ac_table, ncbi_hdr, hdr_len, chr_nb, nb_len);
    if (options.debug) {
      char *nb = hash_table_lookup(ac_table, ncbi_hdr, hdr_len);
      fprintf (stderr, "Hash table value for %s (len = %d) is chr%s (len = %d)\n", ncbi_hdr, hdr_len, nb, nb_len);
    }
  }
  return 0;
}

int 
process_bowtie(FILE *input, char *iFile) 
{
  char *s, *res, *buf;
  size_t bLen = LINE_SIZE;
  unsigned int k = 0;

  if (iFile == NULL) {
    iFile = malloc(6 * sizeof(char)); 
    strcpy(iFile, "stdin");
  }
  if (options.debug) 
    fprintf (stderr, "processing file: %s\n", iFile);
  if ((s = malloc(bLen * sizeof(char))) == NULL) {
    perror("process_bowtie: malloc");
    exit(1);
  }
  while ((res = fgets(s, (int) bLen, input)) != NULL) {
    char sc[SCORE_MAX] = "";
    char strand = '\0';
    char ac[HDR_MAX] = "";
    char start[POS_MAX] = "";
    char tag[TAG_MAX] = "";
    size_t cLen = strlen(s);
    int i = 0;

    k++;
    //fprintf (stderr, "LINE: %d\n", k);
    while (cLen + 1 == bLen && s[cLen - 1] != '\n') {
      bLen *= 2;
      if ((s = realloc(s, bLen)) == NULL) {
        perror("process_bowtie: realloc");
        exit(1);
      }
      res = fgets(s + cLen, (int) (bLen - cLen), input);
      cLen = strlen(s);
    }
    if (s[cLen - 1] == '\n')
      s[cLen - 1] = 0;

    buf = s;
    /* Get bowtie fields */
   /* TAG Score/Name */
    i = 0;
    while (*buf != 0 && !isspace(*buf)) {
      if (i >= SCORE_MAX) {
        fprintf(stderr, "Match score too long \"%s\" \n", buf);
        fclose(input);
        exit(1);
      }
      sc[i++] = *buf++;
    }
    if (i < SCORE_MAX)
      sc[i] = 0;
    while (isspace(*buf))
      buf++;
    /* Strand */
    strand = *buf++;
    while (isspace(*buf))
      buf++;
   /* Chromosome AC Header */
    i = 0;
    while (*buf != 0 && !isspace(*buf)) {
      if (i >= HDR_MAX) {
        fprintf(stderr, "AC too long \"%s\" \n", buf);
        fclose(input);
        exit(1);
      }
      ac[i++] = *buf++;
    }
    if (i < HDR_MAX)
      ac[i] = 0;
    while (isspace(*buf))
      buf++;
   /* Start Position */
    i = 0;
    while (isdigit(*buf)) {
      if (i >= POS_MAX) {
        fprintf(stderr, "Start pos too long \"%s\" \n", buf);
        exit(1);
      }
      start[i++] = *buf++;
    }
    if (i < POS_MAX)
      start[i] = 0;
    while (isspace(*buf))
      buf++;
    /* TAG */
    i = 0;
    while (*buf != 0 && !isspace(*buf)) {
      if (i >= TAG_MAX) {
        fprintf(stderr, "Tag too long \"%s\" \n", buf);
        fclose(input);
        exit(1);
      }
      tag[i++] = *buf++;
    }
    if (i < TAG_MAX)
      tag[i] = 0;
/*
    while (isspace(*buf))
      buf++;
*/
#ifdef DEBUG   
    printf("%s\t%c\t%s\t%s\t%s\n", sc, strand, ac, start, tag );
#endif
    if (strand == '-') {
      reverse(tag);
      complement(tag);
    }
    int pos = atoi(start);
    int hdr_len = strlen(ac) + 1;
    int end = pos + tagLen;
    char *nb = hash_table_lookup(ac_table, ac, hdr_len);
    //fprintf (stderr, "Hash table value for %s is chr%s\n", ac, nb);
    /* Convert to BED format */ 
    if (options.mism) {
      int s = 0;
      s = misMatch;
      printf("chr%s\t%d\t%d\t%s\t%d\t%c\n", nb, pos, end, tag, s, strand); 
    } else {
      if (options.norm) {
        int score = atoi(sc)/options.norm;
        printf("chr%s\t%d\t%d\t%s\t%d\t%c\n", nb, pos, end, tag, score, strand);
      } else {
        //printf("%s\t%d\t%s\t%s\t%s\t%c\n", ac, pos, end, tag, sc, strand);
        printf("chr%s\t%d\t%d\t%s\t%s\t%c\n", nb, pos, end, tag, sc, strand);
      }
    }
  }
  if (options.debug) {
    fprintf (stderr, "Done!\n");
    fprintf(stderr, "Numer of read lines %d\n", k);
  }
  return 0;
}

int main(int argc, char *argv[])
{
#ifdef DEBUG
  mcheck(NULL);
  mtrace();
#endif
  FILE *bowtie;
  options.colSep = "\t";

  while (1) {
    int c = getopt(argc, argv, "s:i:l:t:n:m:dh");
    if (c == -1)
      break;
    switch (c) {
    case 'd':
      options.debug = 1;
      break;
    case 'i':
      options.dbPath = optarg;
      options.db = 1;
      break;
    case 'l':
      tagLen = atoi(optarg);
      break;
    case 'h':
      options.help = 1;
      break;
    case 'm':
      misMatch = atoi(optarg);
      options.mism = 1;
      break;
    case 'n':
      options.norm = atoi(optarg);
      break;
    case 's':
      Species = optarg;
      break;
    default:
      printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }

  if (optind > argc || options.help == 1 || Species == NULL || tagLen == 0) {
    fprintf(stderr, "Usage: %s [options] -s <s_assembly (e.g. hg18)> -l <taglen> [<] <Bowtie file|stdin>\n"
             "      where options are:\n"
             "  \t\t -d             Produce debug information\n"
             "  \t\t -h             Show this help text\n"
             "  \t\t -m <mismatch>  Tag mismatch\n"
             "  \t\t -n             Scaling correction factor for score values\n"
             "  \t\t -i <path>      Use <path> to locate the chr_hdr file (default is /home/local/db/genome)\n"
             "\n\tConvert output of bowtie into BED format.\n"
             "\n\tThe score value is reported in the first field of the Bowtie output file.\n"
             "\tIf the tag mismatch is set (-m option), the score is defined as the tag mismatch.\n\n",
             argv[0]);
    return 1;
  }

  //printf("argc %d optind %d\n", argc, optind);

  if (argc > optind) {
    if(!strcmp(argv[optind],"-")) {
      bowtie = stdin;
    } else {
      bowtie = fopen(argv[optind], "r");
      if (NULL == bowtie) {
        fprintf(stderr, "Unable to open '%s': %s(%d)\n",
                argv[optind], strerror(errno), errno);
        exit(EXIT_FAILURE);
      }
      if (options.debug)
        fprintf(stderr, "Processing file %s\n", argv[optind]);
    }
  } else {
    bowtie = stdin;
  }
  if (options.debug) {
    fprintf(stderr, " Arguments:\n");
    if (options.score)
      fprintf(stderr, "Score file: %s\n", options.scFile); 
    if (options.mism)
      fprintf(stderr, "Tag mismatch: %d\n", misMatch); 
    fprintf(stderr, "Bowtie file: %s\n", argv[optind]); 
    if (options.db)
      fprintf(stderr, "DB path for locating the chr_NC_gi file: %s\n", options.dbPath); 
    fprintf(stderr, "Species assembly: %s\n", Species);
  }

  if (process_ac() == 0) {
    if (options.debug)
      fprintf(stderr, "HASH Table for chromosome access identifier initialized\n");
  } else {
    return 1;
  }
  if (process_bowtie(bowtie, argv[optind++]) != 0) {
    return 1;
  }
  return 0;
}
