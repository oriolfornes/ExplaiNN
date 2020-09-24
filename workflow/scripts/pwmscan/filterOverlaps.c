/*
  filterOverlaps.c

  Filter overlapping matches.
  The program locates and filters out for regions or matches that overlap
  in BED-formatted genomic regions.
  
  # Arguments:
  # BED region length (rLen)

  Giovanna Ambrosini, EPFL/ISREC, Giovanna.Ambrosini@epfl.ch

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
#include <sys/stat.h>
#ifdef DEBUG
#include <mcheck.h>
#endif

/*#define BUF_SIZE 4096 */
#define BUF_SIZE 8192
#define LINE_SIZE 1024
#define SEQ_ID  256
#define POS_MAX 16
#define SCORE_MAX 12
#define TAG_MAX 128

typedef struct _options_t {
  int help;
  int debug;
} options_t;

static options_t options;

typedef struct _bedline_t {
  char seq_id[SEQ_ID];
  unsigned long *start;
  unsigned long *end;
  char **tag;
  int *score;
  char *strand;
  int *pflag;
} bedline_t, *bedline_p_t;

bedline_t bed_reg;

int rLen = 0;

void 
filter_regions(int len)
{
  int max_score = 0;
  int win = 0;
  int i, k;

  for (i = 0; i < len; i++) {
    max_score = bed_reg.score[i];
    if (rLen == 0) {
      win = (int)(bed_reg.end[i] - bed_reg.start[i]);
    } else {
      win = rLen;
    }
    //printf ("i %d : pos %lu score %d (maxscore)\n", i, bed_reg.start[i], max_score);
    for (k = i + 1; (unsigned long)bed_reg.start[k] <= (unsigned long)(bed_reg.start[i] + (unsigned long)win) && k < len; k++) {
    //printf ("FOUND overlapping region k %d : pos %lu score %d (max=%d)\n", k, bed_reg.start[k], bed_reg.score[k], max_score);
      if (bed_reg.score[k] > max_score) {
        bed_reg.pflag[i] = 0;
        max_score = bed_reg.score[k];
        i = k;
      } else {
        bed_reg.pflag[k] = 0;
      }
    }
    i = k - 1;
  }
  /* Print out non-overlapping regions (pflag = 1) */
  for (i = 0; i < len; i++) {
#ifdef DEBUG
    fprintf(stderr,"dbg: %s\t%lu\t%lu\t%s\t%d\t%c\n", bed_reg.seq_id, bed_reg.start[i], bed_reg.end[i], bed_reg.tag[i], bed_reg.score[i], bed_reg.strand[i]); 
#endif
    if (bed_reg.pflag[i]) {
      printf("%s\t%lu\t%lu\t%s\t%d\t%c\n", bed_reg.seq_id, bed_reg.start[i], bed_reg.end[i], bed_reg.tag[i], bed_reg.score[i], bed_reg.strand[i]); 
    }
  }
}

int
process_bed(FILE *input, char *iFile) 
{
  char seq_id_prev[SEQ_ID] = "";
  unsigned long start, end;
  int score;
  size_t mLen = BUF_SIZE;
  char *s, *res, *buf;
  size_t bLen = LINE_SIZE;
  unsigned int k = 0;

  if (options.debug && input != stdin) {
    char sort_cmd[1024] = "sort -s -c -k1,1 -k2,2n ";
    fprintf(stderr, "Check whether file %s is properly sorted\n", iFile);
      if (strcat(sort_cmd, iFile) == NULL) {
        fprintf(stderr, "strcat failed\n");
        return 1;
      }
    if (strcat(sort_cmd, " 2>/tmp/sortcheck.out") == NULL) {
      fprintf(stderr, "strcat failed\n");
      return 1;
    }
    fprintf(stderr, "executing : %s\n", sort_cmd);
    int sys_code = system(sort_cmd);
    if (sys_code != 0) {
      fprintf(stderr, "system command failed\n");
      return 1;
    }
    struct stat file_status;
    if(stat("/tmp/sortcheck.out", &file_status) != 0){
      fprintf(stderr, "could not stat\n");
      return 1;
    }
    if (file_status.st_size != 0) {
      fprintf(stderr, "BED file %s is not properly sorted\n", iFile);
      return 1;
    } else {
      system("/bin/rm /tmp/sortcheck.out");
    }
  }
  if ((bed_reg.start = (unsigned long*)calloc(mLen, sizeof(unsigned long))) == NULL) {
    perror("process_bed: malloc");
    exit(1);
  }
  if ((bed_reg.end = (unsigned long*)calloc(mLen, sizeof(unsigned long))) == NULL) {
    perror("process_bed: malloc");
    exit(1);
  }
  if ((bed_reg.score = (int*)calloc(mLen, sizeof(int))) == NULL) {
    perror("process_bed: malloc");
    exit(1);
  }
  if (( bed_reg.tag = (char**)calloc(mLen, sizeof(*(bed_reg.tag)))) == NULL) {
    perror("process_bed: malloc");
    exit(1);
  }
  if ((bed_reg.strand = (char*)calloc(mLen, sizeof(int))) == NULL) {
    perror("process_bed: malloc");
    exit(1);
  }
  if ((bed_reg.pflag = (int*)calloc(mLen, sizeof(int))) == NULL) {
    perror("process_bed: malloc");
    exit(1);
  }
  if ((s = malloc(bLen * sizeof(char))) == NULL) {
    perror("process_bed: malloc");
    exit(1);
  }
#ifdef DEBUG
  int c = 1; 
#endif
  while ((res = fgets(s, (int) bLen, input)) != NULL) {
    size_t cLen = strlen(s);
    char seq_id[SEQ_ID] = "";
    char tag[TAG_MAX] = ""; 
    char s_pos[POS_MAX] = ""; 
    char e_pos[POS_MAX] = ""; 
    char sc[SCORE_MAX] = ""; 
    char strand = '\0';
    unsigned int i = 0;

    while (cLen + 1 == bLen && s[cLen - 1] != '\n') {
      bLen *= 2;
      if ((s = realloc(s, bLen)) == NULL) {
        perror("process_file: realloc");
        exit(1);
      }
      res = fgets(s + cLen, (int) (bLen - cLen), input);
      cLen = strlen(s);
    }
    if (s[cLen - 1] == '\n')
      s[cLen - 1] = 0;

    buf = s;
    /* Get BED fields */
    /* Chrom NB */
    while (*buf != 0 && !isspace(*buf)) {
      if (i >= SEQ_ID) {
        fprintf(stderr, "Chrom NB too long \"%s\" \n", buf);
        exit(1);
      }
      seq_id[i++] = *buf++;
    }
    while (isspace(*buf))
      buf++;
    /* Start Position */
    i = 0;
    while (isdigit(*buf)) {
      if (i >= POS_MAX) {
        fprintf(stderr, "Start position too large \"%s\" \n", buf);
        exit(1);
      }
      s_pos[i++] = *buf++;
    }
    s_pos[i] = 0;
    start = (unsigned long)atoi(s_pos);
    while (isspace(*buf))
      buf++;
    /* End Position */
    i = 0;
    while (isdigit(*buf)) { 
      if (i >= POS_MAX) {
        fprintf(stderr, "End position too large \"%s\" \n", buf);
        exit(1);
      }
      e_pos[i++] = *buf++;
    }
    e_pos[i] = 0;
    end = (unsigned long)atoi(e_pos);
    while (isspace(*buf))
      buf++;
    /* Tag */
    i = 0;
    while (*buf != 0 && !isspace(*buf)) {
      if (i >= TAG_MAX) {
        fprintf(stderr, "Tag too long \"%s\" \n", buf);
        fclose(input);
        exit(1);
      }
      tag[i++] = *buf++;
    }
    tag[i] = 0;
    while (isspace(*buf))
      buf++;
    /* Score */
    i = 0;
    while (isdigit(*buf) || *buf == '+' || *buf == '-') { 
      if (i >= SCORE_MAX) {
        fprintf(stderr, "Score too large \"%s\" \n", buf);
        exit(1);
      }
      sc[i++] = *buf++;
    }
    sc[i] = 0;
    score = atoi(sc);
    while (isspace(*buf))
      buf++;
    /* Strand */
    strand = *buf++;
    while (isspace(*buf))
      buf++;

#ifdef DEBUG
    printf(" [%d] Chr nb: %s   Start: %lu  End: %lu  Tag: %s  Score: %d Strand: %c\n", c++, seq_id, start, end, tag, score, strand);
#endif
    if (k >= mLen - 1) {
      mLen *= 2;
#ifdef DEBUG
      fprintf(stderr, "reallocating memory for bed_reg.start bed_reg.end bed_reg.tag bed_reg.score (k=%d, size=%d)\n", k, (int)mLen);
#endif
      if ((bed_reg.start = (unsigned long *)realloc(bed_reg.start, mLen * sizeof(unsigned long))) == NULL) {
    	perror("process_bed: realloc");
    	exit(1);
      }
      if ((bed_reg.end = (unsigned long *)realloc(bed_reg.end, mLen * sizeof(unsigned long))) == NULL) {
    	perror("process_bed: realloc");
    	exit(1);
      }
      if ((bed_reg.score = (int *)realloc(bed_reg.score, mLen * sizeof(int))) == NULL) {
    	perror("process_bed: realloc");
    	exit(1);
      }
      if ((bed_reg.tag = (char**)realloc(bed_reg.tag, mLen * sizeof(*(bed_reg.tag)))) == NULL) {
        perror("process_bed: malloc");
        exit(1);
      }
      if ((bed_reg.strand = (char *)realloc(bed_reg.strand, mLen * sizeof(int))) == NULL) {
    	perror("process_bed: realloc");
    	exit(1);
      }
      if ((bed_reg.pflag = (int *)realloc(bed_reg.pflag, mLen * sizeof(int))) == NULL) {
    	perror("process_bed: realloc");
    	exit(1);
      }
    }
    /* Check Chromosome/Sequence BEGINNING, process previous chromosomal region and printout results*/
    if (strcmp(seq_id, seq_id_prev) != 0) {
      filter_regions((int)k); 
      strcpy(seq_id_prev, seq_id);
      k = 0;
    }
    strcpy(bed_reg.seq_id, seq_id);
    bed_reg.tag[k] = malloc(strlen(tag) + 1);
    strcpy(bed_reg.tag[k], tag);
    bed_reg.start[k] = start;
    bed_reg.end[k] = end;
    bed_reg.score[k] = score;
    bed_reg.strand[k] = strand;
    bed_reg.pflag[k] = 1;
    k++;
  } /* End of While */
  /* Filter overlapping regions for last chromosome/sequence */ 
  filter_regions((int)k); 
  if (input != stdin) {
    fclose(input);
  }
  return 0;
}

int
main(int argc, char *argv[])
{
#ifdef DEBUG
  mcheck(NULL);
  mtrace();
#endif
  FILE *input;

  while (1) {
    int c = getopt(argc, argv, "dhl:v:");
    if (c == -1)
      break;
    switch (c) {
      case 'd':
        options.debug = 1;
        break;
      case 'h':
        options.help = 1;
        break;
      case 'l':
        rLen = atoi(optarg);
	break;
      default:
        printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }
  if (optind > argc || options.help == 1) {
    fprintf(stderr, "Usage: %s [options] [-l <len>] [<] <BED file>\n"
             "      where options are:\n"
	     "  \t\t -d     Produce debug information and check BED file\n"
	     "  \t\t -h     Show this help text\n"
	     "  \t\t -l     BED Region length (default is %d)\n"
	     "\n\tFilters out overlapping matches or regions represented in BED or BED-like format.\n"
	     "\n\tIf regions are of fixed size, their length must be set via the -l <len> option.\n"
	     "\tThe BED input file MUST BE sorted by sequence name (or chromosome id), position, and strand.\n"
	     "\tOne should check the input BED file with the following command:\n"
	     "\tsort -s -c -k1,1 -k2,2n -k6,6 <BED file>.\n\n"
	     "\tIn debug mode (-d), the program performs the sorting order check.\n\n"
             "\tThe output is a BED-formatted list of non-overlapping matches.\n\n",
	     argv[0], rLen);
      return 1;
  }
  if (argc > optind) {
      if(!strcmp(argv[optind],"-")) {
          input = stdin;
      } else {
          input = fopen(argv[optind], "r");
          if (NULL == input) {
              fprintf(stderr, "Unable to open '%s': %s(%d)\n",
                  argv[optind], strerror(errno), errno);
             exit(EXIT_FAILURE);
          }
          if (options.debug)
             fprintf(stderr, "Processing file %s\n", argv[optind]);
      }
  } else {
      input = stdin;
  }
  if (options.debug) {
    fprintf(stderr, " Arguments:\n");
    if (rLen) {
      fprintf(stderr, " Size of BED regions : %d\n\n", rLen);
    } else {
      fprintf(stderr, " Warning: regions don't have a fixed length (Len = %d)\n\n", rLen);
    }
  }
  if (process_bed(input, argv[optind++]) != 0) {
    return 1;
  }
  return 0;
}
