/*

  Compute the cumulative score distribution of a Position Weight Matrix (PWM) 

  PWMs elements are integer numbers that are calculated as log-likelihoods.

  - If the e-value threshold is set:
            the corresponding score and cut-off percentage are computed.
  - If the cut-off percentage is set:
            the corresponding score and e-value are computed.
  - If the cut-off score is set:
            the corresponding cut-off percentage and e-value are computed.
  
  Giovanna Ambrosini, EPFL/SV, giovanna.ambrosini@epfl.ch

  Copyright (c) 2014
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
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>
#include <ctype.h>
#include <getopt.h>
#include <assert.h>
#include <float.h>
#ifdef DEBUG
#include <mcheck.h>
#endif

#define BUF_SIZE 3072
#define NUCL  4
#define LMAX  100
#define HDR_MAX 132
#define LINE_SIZE 1024
#define MVAL_MAX 16

typedef struct _options_t {
  int help;
  int debug;
  int raw_score;
  int p_value;
  int perc_score;
} options_t;

static options_t options;

static float bg[] = {0.25,0.25,0.25,0.25};

FILE *pwm_in;

int **pwm;
int pwmLen = 10;

float Pvalue;
float Percentage;
int Score;

static int 
read_pwm(FILE *input, char *iFile)
{
  FILE *f = input;
  int l = 0;
  char *s, *res, *buf;
  size_t bLen = LINE_SIZE;
  int p_len = pwmLen;
  char mval[MVAL_MAX] = "";
  int i;

  if (f == NULL) {
    fprintf(stderr, "Could not open file %s: %s(%d)\n",
            iFile, strerror(errno), errno);
    return -1;
  }
  if (options.debug != 0)
    fprintf(stderr, "read_pwm : processing file %s\n", iFile);

  if ((s = malloc(bLen * sizeof(char))) == NULL) {
    perror("process_sga: malloc");
    return(-1);
  }
  /* Read Matrix file line by line */
  while ((res = fgets(s, (int) bLen, f)) != NULL) {
    size_t cLen = strlen(s);

    while (cLen + 1 == bLen && s[cLen - 1] != '\n') {
      bLen *= 2;
      if ((s = realloc(s, bLen)) == NULL) {
        perror("process_file: realloc");
        return(-1);
      }
      res = fgets(s + cLen, (int) (bLen - cLen), f);
      cLen = strlen(s);
    }
    if (s[cLen - 1] == '\n')
      s[cLen - 1] = 0;

    if (s[cLen - 2] == '\r')
      s[cLen - 2] = 0;

    buf = s;
    /* Get PWM fields */
    /* Get first character: if # or > skip line */
    if (*buf == '#' || *buf == '>')
      continue;

   /* Read First column value */
    while (isspace(*buf))
      buf++;
    i = 0;
    while (isdigit(*buf) || *buf == '-') {
      if (i >= MVAL_MAX) {
        fprintf(stderr, "Matrix value is too large \"%s\" \n", buf);
        return(-1);
     }
      mval[i++] = *buf++;
    }
    mval[i] = 0;
    pwm[l][0] = atoi(mval);
    while (isspace(*buf))
      buf++;
    /* Read Second column value */
    i = 0;
    while (isdigit(*buf) || *buf == '-') {
      if (i >= MVAL_MAX) {
        fprintf(stderr, "Matrix value is too large \"%s\" \n", buf);
        return(-1);
     }
      mval[i++] = *buf++;
    }
    mval[i] = 0;
    pwm[l][1] = atoi(mval);
    while (isspace(*buf))
      buf++;
    /* Read Third column value */
    i = 0;
    while (isdigit(*buf) || *buf == '-') {
      if (i >= MVAL_MAX) {
        fprintf(stderr, "Matrix value is too large \"%s\" \n", buf);
        return(-1);
     }
      mval[i++] = *buf++;
    }
    mval[i] = 0;
    pwm[l][2] = atoi(mval);
    while (isspace(*buf))
      buf++;
    /* Read fourth column value */
    i = 0;
    while (isdigit(*buf) || *buf == '-') {
      if (i >= MVAL_MAX) {
        fprintf(stderr, "Matrix value is too large \"%s\" \n", buf);
        return(-1);
     }
     mval[i++] = *buf++;
    }
    mval[i] = 0;
    pwm[l][3] = atoi(mval);
#ifdef DEBUG
    fprintf(stderr, "%3d   %7d   %7d   %7d   %7d\n", l, pwm[l][0], pwm[l][1], pwm[l][2], pwm[l][3]);
#endif
    if (l == p_len-1) {
      /* Reallocate Matrix rows */
      pwm = realloc(pwm, p_len*2*sizeof(int *));
      if (pwm == NULL) {
        fprintf(stderr, "Out of memory\n");
        return 1;
      }
      /* Allocate columns       */
      for (int i = p_len; i < p_len*2; i++) {
        pwm[i] = calloc((size_t)NUCL, sizeof(int));
        if (pwm[i] == NULL) {
          fprintf(stderr, "Out of memory\n");
          return 1;
        }
      }
      p_len *= 2;
    }
    l++;
  }
#ifdef DEBUG
  fprintf(stderr, "PWM length: %d\n", l);
#endif
  fclose(f);
  return l;
}

int
find_max(int *a, int n)
{
  int i, max;

  max = a[0];
  for (i = 1; i < n; i++)
    if (a[i] > max)
      max = a[i];
  return max;
}

int
find_min(int *a, int n)
{
  int i, min;

  min = a[0];
  for (i = 1; i < n; i++)
    if (a[i] < min)
      min = a[i];
  return min;
}

int
max_score(int k)
{
  /* Compute max score of pwm column k */
  int scores[NUCL] = {0};
  int i;
  int max = 0;

  for (i = 0; i < NUCL; i++) {
    scores[i] = pwm[k][i];
  }
  max = find_max(scores, NUCL);
  if (options.debug)
    fprintf(stderr, "max: %d\n", max);
  return max;
}

int
min_score(int k)
{
  /* Compute min score of pwm column k */
  int scores[NUCL] = {0};
  int i;
  int min = 0;

  for (i = 0; i < NUCL; i++) {
    scores[i] = pwm[k][i];
  }
  min = find_min(scores, NUCL);
  if (options.debug)
    fprintf(stderr, "min: %d\n", min);
  return min;
}

static int
process_pwm()
{
  int i;  /* Nucleoitide code  */
  int k;  /* PWM row           */
  int j;  /* PWM Score array   */ 
  int max = 0;
  int tmp_max = 0;
  int min = 0;
  int offset = 0;
  double *p;
  double *q;
  double prob; /* Cumulative Probabiliy */
  double prob_prev;
  double perc;
  double perc_prev;
  int rscore;
  int rscore_prev;

  /* Rescale PWM to set min=0 for each pwm position/row */
  /* Save nex max score and offset */
  for (k = 0; k < pwmLen; k++) {
    min = min_score(k); 
    tmp_max = max_score(k);
    for (i = 0; i < NUCL; i++)
      pwm[k][i] -= min;
    max += tmp_max - min;
    offset -= min;
  }
  if (options.debug) {
    fprintf(stderr, " Rescaled Weight Matrix: \n\n");
    for (k = 0; k < pwmLen; k++) {
      for ( int i = 0; i < NUCL; i++) {
        int mval = pwm[k][i];
        fprintf(stderr, " %7d ", mval);
      }
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
  }
  /* Score range/space goes now from 0 to max */
  /* Score/Probability array                  */
  // printf("PWM MAX = %d\n", max);
  p = (double *)calloc(max + 1, sizeof(double)); 
  if (p == NULL) {
    fprintf(stderr, "Out of memory\n");
    return 1;
  }
  /* Temp Score/Probability array             */
  q = (double *)calloc(max + 1, sizeof(double)); 
  if (q == NULL) {
    fprintf(stderr, "Out of memory\n");
    return 1;
  }
  for (j = 0; j <= max;  j++) {
    p[j] = 0.0;
    q[j] = 0.0;
  }
  max = 0;
  /* Treat first row (first position)      */
  for (i = 0; i < NUCL; i++) { 
    p[pwm[0][i]] += bg[i];
    // printf ("p[%d] = %f\n", pwm[0][i], p[pwm[0][i]]);
  }
  /* max score of first row/position       */
  max += max_score(0);
  // printf("PWM MAX[0] = %d\n", max);
  /* Loop on sub-sequent rows/positions    */
  for (k = 1; k < pwmLen; k++) {
    for (j = 0; j <= max;  j++) {
      if (p[j] != 0) {    /* Only consider PWM scores     */
        for (i = 0; i < NUCL; i++) {
            /* Update PWM score for position k and base i */
            // printf ("Update prob p[%d] = %.2e    %15.14f\n", j, p[j], p[j]);
            q[j + pwm[k][i]] += p[j]*bg[i];
            // printf ("q[%d] = %.2e\n", j + pwm[k][i], q[j+pwm[k][i]]);
        }
      }
    } 
    max +=  max_score(k); /* Update score range adding max score of row k  */
    // printf("PWM MAX[%d] = %d\n", k, max);
    for (j = 0; j <= max;  j++) {
      p[j] = q[j]; /* Update probabilities for position k                  */
      q[j] = 0;
    }
  }
  /* printf("PWM MAX (final) = %d\n", max); */
  /* Print Cumulative Probability Distribution starting from PWM max score */
  prob = 0;
  prob_prev = 0;
  perc_prev = 100;
  rscore_prev = max - offset;
  for (j = max; j >= 0; j--) {
    if (p[j] != 0) {    /* Only consider PWM scores     */
      prob += p[j];
      rscore = j - offset;
      perc = (double)j/(double)max * 100; 
      //printf ("PERC: %6.2f\n", perc);
      if (options.p_value) {
        if (prob > Pvalue) {
          if ((prob - Pvalue) > (Pvalue - prob_prev)) {
            printf ("SCORE : %6i\tPERC : %6.2f%%\n", rscore_prev, perc_prev);
          } else {
            printf ("SCORE : %6i\tPERC : %6.2f%%\n", rscore, perc);
          }
          return 0;
        }
      } else if (options.raw_score) {
        if (rscore <= Score) {
          printf ("PERC : %6.2f%%\tPVAL : %.2e\n", perc, prob);
          return 0;
        }
      } else if (options.perc_score) {
        if (perc <= Percentage) {
          printf ("SCORE : %6i\tPVAL : %.2e\n", rscore, prob);
          return 0;
        }
      } else {
        printf ("%6i %.2e %6.2f%%\n", rscore, prob, perc);
      }
      prob_prev = prob;
      rscore_prev = rscore;
      perc_prev = perc;
    }
  } 
  return 0;
}

char** str_split(char* a_str, const char a_delim)
{
    char** result = 0;
    size_t count = 0;
    char* tmp = a_str;
    char* last_comma = 0;
    char delim[2];
    delim[0] = a_delim;
    delim[1] = 0;

    /* Count how many elements will be extracted.     */
    while (*tmp) {
        if (a_delim == *tmp) {
            count++;
            last_comma = tmp;
        }
        tmp++;
    }
    /* Add space for trailing token.                  */
    count += last_comma < (a_str + strlen(a_str) - 1);
    /* Add space for terminating null string so caller
       knows where the list of returned strings ends. */
    count++;
    result = malloc(sizeof(char*) *count);

    if (result) {
        size_t idx  = 0;
        char* token = strtok(a_str, delim);

        while (token) {
            assert(idx < count);
            *(result + idx++) = strdup(token);
            token = strtok(0, delim);
        }
        assert(idx == count - 1);
        *(result + idx) = 0;
    }
    return result;
}

int
main(int argc, char *argv[])
{
  /*char *pwmFile = NULL; */
  char *bgProb = NULL;
  char** tokens;
  int i = 0;
  int j = 0;

  int option_index = 0;
  static struct option long_options[] =
      {
          /* These options may or may not set a flag.
             We distinguish them by their indices. */
          {"debug",   no_argument,       0, 'd'},
          {"help",    no_argument,       0, 'h'},
          {"bg",      required_argument, 0, 'b'},
          {"eval",    required_argument, 0, 'e'},
          {"perc",    required_argument, 0, 'p'},
          {"score",   required_argument, 0, 's'},
          {0, 0, 0, 0}
      };
  
#ifdef DEBUG
  mcheck(NULL);
  mtrace();
#endif
  while (1) {
    int c = getopt_long(argc, argv, "dhb:e:p:s:", long_options, &option_index);
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
      bgProb = optarg;
      break;
    case 'e':
      Pvalue = atof(optarg);
      options.p_value = 1;
      break;
    case 'p':
      Percentage = atof(optarg);
      options.perc_score = 1;
      break;
    case 's':
      Score = atoi(optarg);
      options.raw_score = 1;
      break;
    case '?':
      break;
    default:
      printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }
  /* printf ("optind: %d  argc:%d\n", optind, argc);  */
  if (optind > argc || options.help) {
    fprintf(stderr,
	    "Usage: %s [options] [<] <pwm_file>\n"
	    "   where options are:\n"
	    "     -d[--debug]             Produce debugging output\n"
	    "     -h[--help]              Show this stuff\n"
	    "     -b[--bg] <bg freq>      Set the background nucleotide frequencies <bg freq>: 0.25,0.25,0.25,0.25\n"
	    "                             Note that nucleotide frequencies (<bg freq>) MUST BE comma-separated.\n"
	    "     -e[--eval]  <p-value>   Compute raw score and percentage cut-offs corresponding to the given <p-value>\n"
	    "     -p[--perc]  <perc co>   Compute raw score and p-value cut-offs corresponding to the given <perc co>\n"
	    "     -s[--score] <score>     Compute p-value and percentage cut-offs corresponding to the given <score>\n"
	    "\n\tCompute the cumulative score distribution of an integer position weight matrix (<pwm_file>) or PWM.\n"
	    "\tThe PWM weights are integer numbers calculated as log likelihoods (or log-odds).\n"
	    "\tIf the p-value threshold is set, the corresponding score and percentage cut-off values are computed.\n"
	    "\tIf the cut-off percentage is set, the corresponding score and p-value cut-off values are computed.\n"
	    "\tIf the cut-off score is set, the corresponding p-value and percentage cut-off values are computed.\n\n",
	    argv[0]);
    return 1;
  }
  if (argc > optind) {
      if(!strcmp(argv[optind],"-")) {
          pwm_in = stdin;
      } else {
          pwm_in = fopen(argv[optind], "r");
          if (pwm_in == NULL) {
              fprintf(stderr, "Unable to open '%s': %s(%d)\n",
                  argv[optind], strerror(errno), errno);
             exit(EXIT_FAILURE);
          }
          if (options.debug)
             fprintf(stderr, "Processing file %s\n", argv[optind]);
      }
  } else {
      pwm_in = stdin;
  }

  /* Allocate initial memory for the PWM       */
  pwm = (int **)calloc((size_t)pwmLen, sizeof(int *));   /* Allocate rows (PWM length) */
  if (pwm == NULL) {
    fprintf(stderr, "Could not allocate matrix array: %s(%d)\n",
        strerror(errno), errno);
    return 1;
  }
  for (i = 0; i < pwmLen; i++) {
    pwm[i] = calloc((size_t)NUCL, sizeof(int));          /* Allocate columns (NUCL=4)  */
    if (pwm[i] == NULL) {
      fprintf(stderr, "Out of memory\n");
      return 1;
    }
  }

  /* Read Matrix from file  and compute pwmLen */
  if ((pwmLen = read_pwm(pwm_in, argv[optind])) <= 0)
    return 1;

  /* Treat background nucleotide frequencies   */
  if (bgProb != NULL) {
    tokens = str_split(bgProb, ',');
    if (tokens) {
      int i;
      for (i = 0; *(tokens + i); i++) {
        bg[i] =  atof(*(tokens + i));
        free(*(tokens + i));
      }
      if (i != 4) {
        fprintf(stderr, "Number of TOKENS: %d\n", i);
        fprintf(stderr, "Please, specify correct library-dependent nucleotide frequencies <bg freq>: they MUST BE comma-separated!\n");
        exit(1);
      }
      free(tokens);
    }
  }
  if (options.debug != 0) {
    if (pwm_in != stdin) {
      fprintf(stderr, "PWM File : %s\n", argv[optind]);
    } 
    fprintf(stderr, "PWM length: %d\n", pwmLen);
    fprintf(stderr, "Weight Matrix: \n\n");
    for (j = 0; j < pwmLen; j++) {
      for (i = 0; i < NUCL; i++) {
        int mval = pwm[j][i];
        fprintf(stderr, " %7d ", mval);
      }
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
    if (bgProb != NULL) 
      fprintf(stderr, "Background nucleotide frequencies:[%s]\n", bgProb); 
    else
      fprintf(stderr, "Background nucleotide frequencies:\n"); 
    for (i = 0; i < NUCL; i++) {
      fprintf(stderr, "bg[%i] = %f ", i, bg[i]);
    }
    fprintf(stderr, "\n\n");
    if (options.p_value) {
      fprintf(stderr, "P-value cut-off: %f\n", Pvalue);
    }
    if (options.raw_score) {
      fprintf(stderr, "Raw Score cut-off: %d\n", Score);
    }
    if (options.perc_score) {
      fprintf(stderr, "Percentage cut-off: %f\n", Percentage);
    }
    fprintf(stderr, "\n");
  }
  
  if (process_pwm() != 0)
    return 1;
  
  for (i = 0; i < pwmLen; i++)
    free(pwm[i]);
  free(pwm);

  return 0;
}
