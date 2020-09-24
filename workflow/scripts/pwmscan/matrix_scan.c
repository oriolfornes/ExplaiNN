/*
  Scan a DNA sequence file for matches to an INTEGER position weight matrix (PWM)
  
  The DNA sequence is in FASTA format.

  # Arguments:

     # Matrix File
     # Cut-off score (integer)
     # Search mode: both strands/forward [def: both]
     # Word index length 
     # Background model (base composition) 
       a comma-separated list of four numbers, e.g. 25,25,25,25,
       internally normalized to probabilities [def: 0.25,0.25,0.25,0.25] 
     # Sequence File
  
  The matrix format is integer log-odds, where each column represents a 
  nucleotide base in the following order: A, C, G, T. 
  
  The program output a list of PWM matches in BED-like format.

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
#include <sys/types.h>
#include <unistd.h>
#include <ctype.h>
#include <getopt.h>
#include <assert.h>
#include <limits.h>
#ifdef DEBUG
#include <mcheck.h>
#endif

#define BUF_SIZE 4194304 /* 4MB */
#define THIRTY_TWO_MEG 0x2000000ULL
#define LINE_SIZE 1024
#define NUCL  5
#define LMAX  100
#define HDR_MAX 256
#define MVAL_MAX 16

typedef struct _options_t {
  int help;
  int debug;
  int forward;
} options_t;

static options_t options;

static char nucleotide[] = {'N','A','C','G','T'};
static float bgcomp[] = {0.25,0.25,0.25,0.25, 0.0};

typedef struct _seq_t {
  char *hdr;
  short int *seq;
  unsigned int len;
} seq_t, *seq_p_t;

typedef struct _arr_idx_t {
  float value;
  int index;
} arr_idx_t, *arr_idx_p_t;

FILE *fasta_in;

int **pwm;
int **pwm_r;      /* reverse PWM  */
int pwmLen = 10;
int wordLen = 7;

/* Score arrays  */
int *ScoreF;
int *ScoreR;

int cutOff = INT_MIN;
int Offset = 0;

/* Array z is used to compute the next word index (seq[2...j+1])  */
unsigned int *z;

/* PWMs Core Regions (set by define_search_strategy() function)   */
/* In case the PWM is longer than the Word index, we must define  */
/* a core region within the PWM such that it minimizes the sum of */ 
/* weigths for rapid drop-off. The lateral positions are ranked   */
/* by weigth in decreasing order of importance.                   */
int Bfw = 0;     /* Beginning forward Core Region                 */
int Efw = 0;     /* End forward Core Region                       */
int Brv = 0;     /* Beginning reverse Core Region                 */
int Erv = 0;     /* End reverse Core Region                       */

int *Rfw;        /* Ranked index array (for lateral positions) FW */
int *Rrv;        /* Ranked index array (for lateral positions) RV */

/* Number of Pipe delimiters in the FASTA header after which the seq ID starts */
int nbPipes = 2;

/* Input process functions  */
static int 
read_pwm(char *iFile)
{
  FILE *f = fopen(iFile, "r");
  int l = 0;
  char *s, *res, *buf;
  size_t bLen = LINE_SIZE;
  int p_len = pwmLen + 1;
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
    exit(1);
  }
  /* Read Matrix file line by line */
  while ((res = fgets(s, (int) bLen, f)) != NULL) {
    size_t cLen = strlen(s);

    while (cLen + 1 == bLen && s[cLen - 1] != '\n') {
      bLen *= 2;
      if ((s = realloc(s, bLen)) == NULL) {
        perror("process_file: realloc");
        exit(1);
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

    l++;
    if (l == p_len) {
      /* Reallocate Matrix rows */
      pwm = realloc(pwm, p_len*2*sizeof(int *));
      if (pwm == NULL) {
        fprintf(stderr, "Out of memory\n");
        return 1;
      }
      /* Allocate columns       */
      for ( int i = p_len; i <= p_len*2; i++) {
        pwm[i] = calloc((size_t)NUCL, sizeof(int));
        if (pwm[i] == NULL) {
          fprintf(stderr, "Out of memory\n");
          return 1;
        }
      }
      p_len *= 2;
    }
    /* Read First column value */ 
    while (isspace(*buf))
      buf++;
    i = 0;
    while (isdigit(*buf) || *buf == '-') {
      if (i >= MVAL_MAX) {
        fprintf(stderr, "Matrix value is too large \"%s\" \n", buf);
        exit(1);
     }
      mval[i++] = *buf++;
    }
    mval[i] = 0;
    pwm[l][1] = atoi(mval);
    while (isspace(*buf))
      buf++;
    /* Read Second column value */ 
    i = 0;
    while (isdigit(*buf) || *buf == '-') {
      if (i >= MVAL_MAX) {
        fprintf(stderr, "Matrix value is too large \"%s\" \n", buf);
        exit(1);
     }
      mval[i++] = *buf++;
    }
    mval[i] = 0;
    pwm[l][2] = atoi(mval);
    while (isspace(*buf))
      buf++;
    /* Read Third column value */ 
    i = 0;
    while (isdigit(*buf) || *buf == '-') {
      if (i >= MVAL_MAX) {
        fprintf(stderr, "Matrix value is too large \"%s\" \n", buf);
        exit(1);
     }
      mval[i++] = *buf++;
    }
    mval[i] = 0;
    pwm[l][3] = atoi(mval);
    while (isspace(*buf))
      buf++;
    /* Read fourth column value */ 
    i = 0;
    while (isdigit(*buf) || *buf == '-') {
      if (i >= MVAL_MAX) {
        fprintf(stderr, "Matrix value is too large \"%s\" \n", buf);
        exit(1);
     }
      mval[i++] = *buf++;
    }
    mval[i] = 0;
    pwm[l][4] = atoi(mval);
#ifdef DEBUG
    fprintf(stderr, "%3d   %7d   %7d   %7d   %7d\n", l, pwm[l][1], pwm[l][2], pwm[l][3], pwm[l][4]);
#endif
  }
  fclose(f);
#ifdef DEBUG
  fprintf(stderr, "PWM length: %d\n", l);
#endif
  if (!options.forward) {
    p_len = pwmLen + 1;
    /* Make reverse-complement PWM */
    for (int k = 1; k <= l; k++) { 
      for (i = 1; i < NUCL; i++) 
        pwm_r[k][i] = pwm[l-k+1][NUCL-i];

      if (k == p_len -1) {
      /* Reallocate Matrix rows */
        pwm_r = realloc(pwm_r, p_len*2*sizeof(int *));
        if (pwm_r == NULL) {
          fprintf(stderr, "Out of memory\n");
          return 1;
        }
        /* Allocate columns       */
        for ( int i = p_len; i <= p_len*2; i++) {
          pwm_r[i] = calloc((size_t)NUCL, sizeof(int));
          if (pwm_r[i] == NULL) {
            fprintf(stderr, "Out of memory\n");
            return 1;
          }
        }
        p_len *= 2;
      }
    }
  }
  if ( l == 0 ) return -1;
  return l;
}

static int
find_max(int *a, int n)
{
  int i, max;

  max = a[0];
  for (i = 1; i < n; i++)
    if (a[i] > max)
      max = a[i];
  return max;
}

static int
max_score(int **m, int k)
{
  /* Compute max score of pwm column k */
  int scores[NUCL-1] = {0};
  int i;
  int max = 0;

  for (i = 0; i < NUCL-1; i++) {
    scores[i] = m[k][i+1];
  }
  max = find_max(scores, NUCL-1);
  if (options.debug)
    fprintf(stderr, "max: %d\n", max);
  return max;
}

static void
process_pwm() { 
  /* Rescale weights to have zero as a maximum value at each position */
  /* Compute Offset and re-define cutOff                              */ 
  int max;
  for (int k = 1; k <= pwmLen; k++) {
    max = max_score(pwm, k);
    for (int i = 1; i < NUCL; i++)
      pwm[k][i] -= max;
    Offset += max;
  }
  cutOff = cutOff - Offset; 
  if (options.debug)
    fprintf(stderr, "rescaled cutOff: %d\n", cutOff);
  if (!options.forward) {
  /* Rescale reverse PWM */
    for (int k = 1; k <= pwmLen; k++) {
      max = max_score(pwm_r, k);
      for (int i = 1; i < NUCL; i++)
        pwm_r[k][i] -= max;
    }
  }
  if (options.debug) {
    fprintf(stderr, "Re-scaled Weight Matrix: original representation \n\n");
    for (int k = 1; k <= pwmLen; k++) {
      for ( int i = 1; i < NUCL; i++) {
        int mval = pwm[k][i];
        fprintf(stderr, " %7d ", mval);
      }
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
    if (!options.forward) {
      fprintf(stderr, "Re-scaled Reverse Weight Matrix:\n\n");
      for (int k = 1; k <= pwmLen; k++) {
        for ( int i = 1; i < NUCL; i++) {
          int mval = pwm_r[k][i];
          fprintf(stderr, " %7d ", mval);
        }
        fprintf(stderr, "\n");
      }
      fprintf(stderr, "\n");
    }
  }
}

static void
process_bgcomp() {
  float sum = 0;
  for (int i = 0; i < NUCL-1; i++)
    sum += bgcomp[i];
  for (int i = NUCL-1; i > 0; i--) 
    bgcomp[i] = bgcomp[i-1]/sum;
}

static int 
compfunc(const void *e1, const void *e2) {
  arr_idx_p_t first = (arr_idx_p_t) e1;
  arr_idx_p_t second = (arr_idx_p_t) e2;
  if ((*first).value > (*second).value) return  1;
  if ((*first).value < (*second).value) return -1;
  return 0;
}

/* Prepare Word index and Score tables and define search strategy */
static void
define_search_strategy() { 
  /* Definition of a core region within the PWM */
  /* and ranking of the positions outside the   */
  /* the core region by decreasing importance   */
  float wf[pwmLen+1];
  arr_idx_t wfobj[pwmLen+1];
  /* Allocate forward ranked index array */
  Rfw = (int *) calloc((size_t)pwmLen+1, sizeof(int)); 
  wf[0] = 1.0;
  for (int k = 1; k <= pwmLen; k++) {
    wf[k] = 0.0;
    for (int i = 1; i < NUCL; i++)
      wf[k] += bgcomp[i]*pwm[k][i]; 
  }
  /* Determine core region [Bfw-Efw] minimizing the sum of weights       */
  float x = 0.0;
  for (int j = 1; j <= wordLen; j++)
    x += wf[j];
  float min = x;
  int pos = wordLen;
  for (int j = wordLen+1; j <= pwmLen; j++) {
    x = x - wf[j-wordLen] + wf[j];
    if (x < min) {
      min = x;
      pos = j;
    }
  }
  Bfw = pos - wordLen + 1;
  Efw = pos;
  if (options.debug)
    fprintf (stderr, "Core region FW: from %d to %d\n", Bfw, Efw);
  /* Mask core region and rank lateral positions by weigth               */
  for (int j = Bfw; j <= Efw; j++)
    wf[j] = 1.0;
  /* Fill weight array structure  */
  for (int j = 0; j <= pwmLen; j++) {
    wfobj[j].value = wf[j];
    wfobj[j].index = j;
    if (options.debug)
      fprintf(stderr, "%f ", wf[j]);
  }
  if (options.debug)
    fprintf(stderr, "\n");
  qsort(wfobj, pwmLen+1, sizeof(wfobj[0]), compfunc);
  /* Extract sorted index array   */ 
  for (int j = 0; j <= pwmLen; j++) {
    Rfw[j] = wfobj[j].index;
    if (options.debug)
      fprintf (stderr, " %d ", Rfw[j]);
  }
  if (options.debug)
    fprintf(stderr, "\n");
  if (!options.forward) {  /* Reverse PWM  */
    float wr[pwmLen+1]; 
    arr_idx_t wrobj[pwmLen+1];
    /* Allocate reverse ranked index array */
    Rrv = (int *) calloc((size_t)pwmLen+1, sizeof(int)); 
    wr[0] = 1.0;
    for (int k = 1; k <= pwmLen; k++) {
      wr[k] = 0.0;
      for (int i = 1; i < NUCL; i++)
        wr[k] += bgcomp[i]*pwm_r[k][i]; 
    }
    /* Determine core region [Brv-Erv] minimizing the sum of weights       */
    float x = 0.0;
    for (int j = 1; j <= wordLen; j++)
      x += wr[j];
    float min = x;
    int pos = wordLen;
    for (int j = wordLen+1; j <= pwmLen; j++) {
      x = x - wr[j-wordLen] + wr[j];
      if (x < min) {
        min = x;
        pos = j;
      }
    }
    Brv = pos - wordLen + 1;
    Erv = pos;
    if (options.debug)
      fprintf (stderr, "Core region RV: from %d to %d\n", Brv, Erv);
    /* Mask core region and rank lateral positions by weigth               */
    for (int j = Brv; j <= Erv; j++)
      wr[j] = 1.0;
    /* Fill weight array structure  */
    for (int j = 0; j <= pwmLen; j++) {
      wrobj[j].value = wr[j];
      wrobj[j].index = j;
      if (options.debug)
        fprintf(stderr, "%f ", wr[j]);
    }
    if (options.debug)
      fprintf(stderr, "\n");
    qsort(wrobj, pwmLen+1, sizeof(wrobj[0]), compfunc);
    /* Extract sorted index array   */ 
    for (int j = 0; j <= pwmLen; j++) {
      Rrv[j] = wrobj[j].index;
      if (options.debug)
        fprintf (stderr, " %d ", Rrv[j]);
    }
    if (options.debug)
      fprintf(stderr, "\n");
  }
}

unsigned int
power(unsigned int base, unsigned int degree)
{
  unsigned int result = 1;
  unsigned int term = base;
  while (degree) {
    if (degree & 1)
      result *= term;
    term *= term;
    degree = degree >> 1;
  }
  return result;
}

static int
make_tables()
{
  /* Make Word index and score tables                                           */
  /* Words of length wordLen are encoded as integers between 1 and 4^(wordLen), */ 
  /* e.g. for wordLen=4 index(AAAA)=1, and index(TTTT)=256.                     */
  /* The PWM scores for these words in forward and reverse orientation are      */
  /* stored in ScoreR and ScoreF (integer array variables).                     */ 
  unsigned int i; /* word index 1..4^(wordLen) */
  int j;
  int n;
  /* Allocate memory for score arrays  */
  unsigned int wsize = power(4, wordLen);
  /* Allocate forward score and navigation link arrays                          */
  if ( (ScoreF = (int *)malloc((size_t)(wsize+1) * sizeof(int))) == NULL ) {
    perror("ScoreF: malloc");
    exit(1);
  }
  /* Variable z is used to set the word idx e.g.: z={0,64,128,192} for wordLen=4*/
  if ( (z = (unsigned int *) calloc((size_t)NUCL, sizeof(unsigned int))) == NULL ) {
    perror("z: calloc");
    exit(1);
  }
  for (j = 0; j < NUCL-1; j++) {
    z[j+1] = j * (power(4, wordLen-1));  
    if (options.debug)
      fprintf(stderr, "z[%d] = %d\n", j+1, z[j+1]);
  }
  /* Allocate word array s[0..wordLen+1]: it stores the sequence (numerical form)*/ 
  int *s = (int *) calloc((size_t)wordLen+1, sizeof(int));
  if (s == NULL) {
    perror("s: calloc");
    exit(1);
  }
  /* Allocate temp word score array xf[0..wordLen+1]/xr[0..wordLen+1]            */ 
  int *xf = (int *) calloc((size_t)wordLen+1, sizeof(int));
  if (xf == NULL) {
    perror("xf: calloc");
    exit(1);
  }
  int *xr = (int *) calloc((size_t)wordLen+1, sizeof(int));
  if (xr == NULL) {
    perror("xr: calloc");
    exit(1);
  }
  if (wordLen != pwmLen) {
    if (options.forward) {
      s[1] = 0;
      n = 1; /* partial word lenght (1..wordLen)      */ 
      xf[0] = 0;
      i = 1; 
      while (n > 0) {
        /* Loop over the entire word index            */ 
        s[n]++;
        /* Compute word score up to wordlen n         */
        xf[n] = xf[n-1] + pwm[Bfw+n-1][s[n]];
        /* Compute word score for the remaining part  */
        for(j = n + 1; j <= wordLen; j++) {
          n++;
          s[n] = 1;   /* set character to A           */
          xf[n] = xf[n-1] + pwm[Bfw+n-1][1];
        }
        /* Set word score and navigation link         */ 
        ScoreF[i] = xf[n]; /*  n=wordLen              */
#ifdef DEBUG
        fprintf(stderr, "%u  ", i);
        for (int k = 1; k <= wordLen; k++)
          fprintf(stderr, "%d ", s[k]);
        fprintf(stderr, "  %8d\n", ScoreF[i]);
#endif
        i++;
        /* Keep decreasing n by 1 while s[n]=T        */
        while(s[n] == 4)
          n--;
      }
    } else { /* Scan both strands                     */
      s[1] = 0;
      n = 1; /* partial word lenght (1..wordLen)      */ 
      xf[0] = 0;
      xr[0] = 0;
      i = 1; 
      /* Allocate reverse score array                 */
      if ( (ScoreR = (int *)malloc((size_t)(wsize+1) * sizeof(int))) == NULL) {
        perror("ScoreF: malloc");
        exit(1);
      }
      while (n > 0) {
        /* Loop over the entire word index            */ 
        s[n]++;
        /* Compute word score up to wordlen n         */
        xf[n] = xf[n-1] + pwm[Bfw+n-1][s[n]];
        xr[n] = xr[n-1] + pwm_r[Brv+n-1][s[n]];
        /* Compute word score for the remaining part  */
        for(j = n + 1; j <= wordLen; j++) {
          n++;
          s[n] = 1;   /* set character to A           */
          xf[n] = xf[n-1] + pwm[Bfw+n-1][1];
          xr[n] = xr[n-1] + pwm_r[Brv+n-1][1];
        }
        /* Set word scores and navigation link        */ 
        ScoreF[i] = xf[n];  /* n=wordLen              */
        ScoreR[i] = xr[n];  /* n=wordLen              */
#ifdef DEBUG
        fprintf(stderr, "%u  ", i);
        for (int k = 1; k <= wordLen; k++)
          fprintf(stderr, "%d ", s[k]);
        fprintf(stderr, "  %8d  %8d\n", ScoreF[i], ScoreR[i]);
#endif
        i++;
        /* Keep decreasing n by 1 while s[n]=T        */
        while(s[n] == 4)
          n--;
      }
    }
  } else { /* PWM lenght and word lenght are the same */
    if (options.forward) {
      s[1] = 0;
      n = 1; /* partial word lenght (1..wordLen)      */ 
      xf[0] = 0;
      i = 1; 
      while (n > 0) {
        /* Loop over the entire word index            */ 
        s[n]++;
        /* Compute word score up to wordlen n         */
        xf[n] = xf[n-1] + pwm[n][s[n]];
        /* Compute word score for the remaining part  */
        for(j = n + 1; j <= wordLen; j++) {
          n++;
          s[n] = 1;   /* set character to A           */
          xf[n] = xf[n-1] + pwm[n][1];
        }
        /* Set word score                             */ 
        ScoreF[i] = xf[n];  /* n=wordLen              */
#ifdef DEBUG
        fprintf(stderr, "%u  ", i);
        for (int k = 1; k <= wordLen; k++)
          fprintf(stderr, "%d ", s[k]);
        fprintf(stderr, "  %8d\n", ScoreF[i]);
#endif
        i++;
        /* Keep decreasing n by 1 while s[n]=T        */
        while(s[n] == 4)
          n--;
      }
    } else { /* Scan both strands                     */
      s[1] = 0;
      n = 1; /* partial word lenght (1..wordLen)      */ 
      xf[0] = 0;
      xr[0] = 0;
      i = 1; 
      /* Allocate reverse score array                 */
      if ( (ScoreR = (int *)malloc((size_t)(wsize+1) * sizeof(int))) == NULL) {
        perror("ScoreF: malloc");
        exit(1);
      }
      while (n > 0) {
        /* Loop over the entire word index            */ 
        s[n]++;
        /* Compute word score up to wordlen n         */
        xf[n] = xf[n-1] + pwm[n][s[n]];
        xr[n] = xr[n-1] + pwm_r[n][s[n]];
        /* Compute word score for the remaining part  */
        for(j = n + 1; j <= wordLen; j++) {
          n++;
          s[n] = 1;   /* set character to A           */
          xf[n] = xf[n-1] + pwm[n][1];
          xr[n] = xr[n-1] + pwm_r[n][1];
        }
        /* Set word scores                            */ 
        ScoreF[i] = xf[n];  /* n=wordLen              */
        ScoreR[i] = xr[n];  /* n=wordLen              */
#ifdef DEBUG
        fprintf(stderr, "%u  ", i);
        for (int k = 1; k <= wordLen; k++)
          fprintf(stderr, "%d ", s[k]);
        fprintf(stderr, "  %8d  %8d\n", ScoreF[i], ScoreR[i]);
#endif
        i++;
        /* Keep decreasing n by 1 while s[n]=T        */
        while(s[n] == 4)
          n--;
      }
    }
  }
  free(s);
  free(xf);
  free(xr);
  return 0;
}

static unsigned int
word_index(seq_p_t seq, unsigned int j)
{
  /* Compute index for word startind at sequence position j */
  unsigned int k = j + wordLen -1;
  unsigned int index = 1;
  int base = 1;
  while(k >= j) {
    index += (seq->seq[k] -1) * base;
    k--;
    base = base << 2;
  } 
  return index; 
}

/* Scanning functions */
static void
scan_seq_1f(seq_p_t seq) 
{
  if (seq->len >= (unsigned long)pwmLen) { /*      Forward Scanning          */
    unsigned int j = 0;
    seq->seq[0] = 0;
    unsigned int i = 0;

    while (j <= seq->len) {    /* Loop through the entire sequence           */
      if (seq->seq[j] == 0) {  /* Check if beginning of the process          */
        int n = 0;             /* or we found an N                           */
        while ((n < pwmLen) && (j < seq->len)) {
          /* This loop serves to find the end position (j) of the next       */
          /* that doesn's contain any N's (i.e. seq[j]=0)                    */
          /* The loop terminates when either a word without N's is found     */ 
          /* (i.e. n=8) or the end of teh sequence (j=seq->len) is reached   */
          j++;
          if (seq->seq[j] == 0) /* Check whether we have N's and keep        */
            n = 0;              /* incrementing j (n stays at 0 if N)        */
          else 
            n++;                /* increment n if found ACGT base            */
        }
        if (n < pwmLen)   /* It is used to break the scanning loop           */
          break;          /* in case the last characters are N's             */
        /* Compute word index */
        i = word_index(seq, j - pwmLen + 1); 
      } else { /* We are not at the beginning of the scanning process        */ 
               /* Compute next word index (after shifting by 1 bp)           */ 
               /* using the previous word index and the z link array         */
        i = ((i -z[seq->seq[j-pwmLen]] -1)<<2) + seq->seq[j];
      }
      /* Check for match (j points to the end of candidate sequence)         */
      int score = ScoreF[i];
      if (score >= cutOff) {
        score = score + Offset;
        unsigned int pos = j;
        printf("%s\t%u\t%u\t", seq->hdr, pos-pwmLen, pos);
        /* print word */
        unsigned int k = 0;
        for (k = j-pwmLen+1; k <= j; k++)
          printf("%c", nucleotide[seq->seq[k]]);   
        /* print score */
        printf("\t%d\t+\n", score);
      }
      /* Move on to the next position                                        */
      j++;
    } /* Scanning loop                                                       */
  }
}

static void
scan_seq_1(seq_p_t seq) 
{
  if (seq->len >= (unsigned long)pwmLen) { /*    Bidirectional Scanning      */
    unsigned int j = 0;
    seq->seq[0] = 0;
    unsigned int i = 0;
    
    while (j <= seq->len) {    /* Loop through the entire sequence           */
      if (seq->seq[j] == 0) {  /* Check if beginning of the process          */
        int n = 0;             /* or we found an N                           */
        while ((n < pwmLen) && (j < seq->len)) {
          /* This loop serves to find the end position (j) of the next       */
          /* that doesn's contain any N's (i.e. seq[j]=0)                    */
          /* The loop terminates when either a word without N's is found     */ 
          /* (i.e. n=8) or the end of teh sequence (j=seq->len) is reached   */
          j++;
          if (seq->seq[j] == 0) /* Check whether we have N's and keep        */
            n = 0;              /* incrementing j (n stays at 0 if N)        */
          else 
            n++;                /* increment n if found ACGT base            */
        }
        if (n < pwmLen)   /* It is used to break the scanning loop           */
          break;          /* in case the last characters are N's             */
        /* Compute word index */
        i = word_index(seq, j - pwmLen + 1); 
      } else { /* We are not at the beginning of the scanning process        */ 
               /* Compute next word index (after shifting by 1 bp)           */ 
               /* using the previous word index and the z link array         */
        i = ((i -z[seq->seq[j-pwmLen]] -1)<<2) + seq->seq[j];
      }
      /* Check for match (j points to the end of candidate sequence)         */
      /* Score in forward direction                                          */
      int score = ScoreF[i];
      if (score >= cutOff) {
        score = score + Offset;
        unsigned int pos = j;
        printf("%s\t%u\t%u\t", seq->hdr, pos-pwmLen, pos);
        /* print word */
        unsigned int k = 0;
        for (k = j-pwmLen+1; k <= j; k++)
          printf("%c", nucleotide[seq->seq[k]]);   
        /* print score */
        printf("\t%d\t+\n", score);
      }
      /* Score in reverse direction                                          */
      score = ScoreR[i];
      if (score >= cutOff) {
        score = score + Offset;
        unsigned int pos = j;
        printf("%s\t%u\t%u\t", seq->hdr, pos-pwmLen, pos);
        /* print word */
        unsigned int k = 0;
        for (k = j; k > j-pwmLen; k--)
          printf("%c", nucleotide[NUCL-seq->seq[k]]);   
        /* print score */
        printf("\t%d\t-\n", score);
      }
      /* Move on to the next position                                        */
      j++;
    } /* Scanning loop                                                       */
  }
}

static void
scan_seq_2f(seq_p_t seq)  /* Word index length is smaller than pwm length    */ 
{
  if (seq->len >= (unsigned long)pwmLen) { /*      Forward Scanning          */
    int diff = pwmLen - wordLen;
    /* Indexes of most relevant PWM positions relative to the end of the PWM */
    int *Ifw = (int *) calloc((size_t)diff, sizeof(int)); 
    for (int k = 0; k < diff; k++) 
      Ifw[k] = Rfw[k] - pwmLen;
    /* Re-define forward core region relative to the end of the PWM          */
    int Bfw_rel = Bfw - pwmLen;
    int Efw_rel = Efw - pwmLen;

    unsigned int j = 0;
    seq->seq[0] = 0;
    unsigned int i = 0;

    while (j <= seq->len) {    /* Loop through the entire sequence           */
      if (seq->seq[j] == 0) {  /* Check if beginning of the process          */
        int n = 0;             /* or we found an N                           */
        while ((n < pwmLen) && (j < seq->len)) {
          /* This loop serves to find the end position (j) of the next       */
          /* that doesn's contain any N's (i.e. seq[j]=0)                    */
          /* The loop terminates when either a word without N's is found     */ 
          /* (i.e. n=8) or the end of teh sequence (j=seq->len) is reached   */
          j++;
          if (seq->seq[j] == 0) /* Check whether we have N's and keep        */
            n = 0;              /* incrementing j (n stays at 0 if N)        */
          else 
            n++;                /* increment n if found ACGT base            */
        }
        if (n < pwmLen)   /* It is used to break the scanning loop           */
          break;          /* in case the last characters are N's             */
        /* Compute word index */
        i = word_index(seq, j + Bfw_rel); 
      } else { /* We are not at the beginning of the scanning process        */ 
               /* Compute next word index (after shifting by 1 bp)           */ 
               /* using the previous word index and the z link array         */
        i = ((i -z[seq->seq[j+Bfw_rel-1]] -1)<<2) + seq->seq[j+Efw_rel];
      }
      /* Check for match (j points to the end of candidate sequence)         */
      int score = ScoreF[i];
      /* Complete score computation with the remaining PWM positions         */
      int k = 0;
      while (score >= cutOff && k < diff) {
        score += pwm[Rfw[k]][seq->seq[j+Ifw[k]]];
        k++;
      }
      if (score >= cutOff) {
        score = score + Offset;
        unsigned int pos = j;
        printf("%s\t%u\t%u\t", seq->hdr, pos-pwmLen, pos);
        /* print word */
        unsigned int k = 0;
        for (k = j-pwmLen+1; k <= j; k++)
          printf("%c", nucleotide[seq->seq[k]]);   
        /* print score */
        printf("\t%d\t+\n", score);
      }
      /* Move on to the next position                                        */
      j++;
    } /* Scanning loop                                                       */
  }
}

static void
scan_seq_2(seq_p_t seq)   /* Word index length is smaller than pwm length    */ 
{
  if (seq->len >= (unsigned long)pwmLen) { /*   Bidirectional Scanning       */
    int diff = pwmLen - wordLen;
    /* Indexes of most relevant PWM positions relative to the end of the PWM */
    int *Ifw = (int *) calloc((size_t)diff, sizeof(int)); 
    int *Irv = (int *) calloc((size_t)diff, sizeof(int)); 
    for (int k = 0; k < diff; k++) { 
      Ifw[k] = Rfw[k] - pwmLen;
      Irv[k] = Rrv[k] - pwmLen;
    }
    /* Re-define forward/rev core regions relative to the end of the PWM     */
    int Bfw_rel = Bfw - pwmLen;
    int Efw_rel = Efw - pwmLen;
    int Brv_rel = Brv - pwmLen;
    int Erv_rel = Erv - pwmLen;

    unsigned int j = 0;
    seq->seq[0] = 0;
    unsigned int ifw = 0;
    unsigned int irv = 0;

    while (j <= seq->len) {    /* Loop through the entire sequence           */
      if (seq->seq[j] == 0) {  /* Check if beginning of the process          */
        int n = 0;             /* or we found an N                           */
        while ((n < pwmLen) && (j < seq->len)) {
          /* This loop serves to find the end position (j) of the next       */
          /* that doesn's contain any N's (i.e. seq[j]=0)                    */
          /* The loop terminates when either a word without N's is found     */ 
          /* (i.e. n=8) or the end of teh sequence (j=seq->len) is reached   */
          j++;
          if (seq->seq[j] == 0) /* Check whether we have N's and keep        */
            n = 0;              /* incrementing j (n stays at 0 if N)        */
          else 
            n++;                /* increment n if found ACGT base            */
        }
        if (n < pwmLen)   /* It is used to break the scanning loop           */
          break;          /* in case the last characters are N's             */
        /* Compute word index */
        ifw = word_index(seq, j + Bfw_rel);
        irv = word_index(seq, j + Brv_rel);
      } else { /* We are not at the beginning of the scanning process        */ 
               /* Compute next word index (after shifting by 1 bp)           */ 
               /* using the previous word index and the z link array         */
        ifw = ((ifw -z[seq->seq[j+Bfw_rel-1]] -1)<<2) + seq->seq[j+Efw_rel];
        irv = ((irv -z[seq->seq[j+Brv_rel-1]] -1)<<2) + seq->seq[j+Erv_rel];
      }
      /* Check for match (j points to the end of candidate sequence)         */

      /* Score in forward direction                                          */
      int score = ScoreF[ifw];
      /* Complete score computation with the remaining PWM positions         */
      int k = 0;
      while (score >= cutOff && k < diff) {
        score += pwm[Rfw[k]][seq->seq[j+Ifw[k]]];
        k++;
      }
      if (score >= cutOff) {
        score = score + Offset;
        unsigned int pos = j;
        printf("%s\t%u\t%u\t", seq->hdr, pos-pwmLen, pos);
        /* print word */
        unsigned int k = 0;
        for (k = j-pwmLen+1; k <= j; k++)
          printf("%c", nucleotide[seq->seq[k]]);   
        /* print score */
        printf("\t%d\t+\n", score);
      }

      /* Score in reverse direction                                          */
      score = ScoreR[irv];
      /* Complete score computation with the remaining PWM positions         */
      k = 0;
      while (score >= cutOff && k < diff) {
        score += pwm_r[Rrv[k]][seq->seq[j+Irv[k]]];
        k++;
      }
      if (score >= cutOff) {
        score = score + Offset;
        unsigned int pos = j;
        printf("%s\t%u\t%u\t", seq->hdr, pos-pwmLen, pos);
        /* print word */
        unsigned int k = 0;
        for (k = j; k > j-pwmLen; k--)
          printf("%c", nucleotide[NUCL-seq->seq[k]]);   
        /* print score */
        printf("\t%d\t-\n", score);
      }
      /* Move on to the next position                                        */
      j++;
    } /* Scanning loop                                                       */
  }
}

/* String parser function */
char** str_split(char* a_str, const char a_delim)
{
    char** result = 0;
    size_t count = 0;
    char* tmp = a_str;
    char* last_comma = 0;
    char delim[2];
    delim[0] = a_delim;
    delim[1] = 0;

    /* Count how many elements will be extracted. */
    while (*tmp) {
        if (a_delim == *tmp) {
            count++;
            last_comma = tmp;
        }
        tmp++;
    }
    /* Add space for trailing token. */
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

/* Process Sequence file - Main Loop */
static int
process_seq(FILE *input, char *iFile)
{
  char buf[BUF_SIZE], *res;
  seq_t seq;
  unsigned int mLen;

  if (input == NULL) {
    FILE *f = fopen(iFile, "r");
    if (f == NULL) {
      fprintf(stderr, "Could not open file %s: %s(%d)\n",
  	    iFile, strerror(errno), errno);
      return -1;
    }
    input = f;
  }
  if (options.debug != 0) {
    if (iFile == NULL)
      fprintf(stderr, "Processing file from STDIN\n");
    else
      fprintf(stderr, "Processing file %s\n", iFile);
  }
  while ((res = fgets(buf, BUF_SIZE, input)) != NULL
	 && buf[0] != '>')
    ;
  if (res == NULL || buf[0] != '>') {
    fprintf(stderr, "Could not find a sequence in file %s\n", iFile);
    if (input != stdin) {
      fclose(input);
    }
    return -1;
  }
  seq.hdr = malloc(HDR_MAX * sizeof(char));
  seq.seq = malloc(THIRTY_TWO_MEG * sizeof(short int));
  mLen = THIRTY_TWO_MEG;
  while (res != NULL) {
    /* Get the header */
    if (buf[0] != '>') {
      fprintf(stderr, "Could not find a sequence header in file %s\n", iFile);
      if (input != stdin) {
        fclose(input);
      }
      return -1;
    }
    char *s = buf;
    s += 1;
    int i = 0;
    while (*s && !isspace(*s)) {
      if (i >= HDR_MAX) {
        fprintf(stderr, "Fasta Header too long \"%s\" in file %s\n", res, iFile);
        fclose(input);
        return -1;
      }
      seq.hdr[i++] = *s++;
    } 
    if (i < HDR_MAX)
      seq.hdr[i] = 0;
    /* Extract sequence identifier from FASTA header               */
    /* The header is expected to have pipe delimiters              */ 
    /* By default, the seq identifier should start after the 2nd   */
    /* pipe (nbPipes=2), and be followed by a space                */ 
    /* The parameter nbPipes can be change by the user via the     */
    /* command line arguments                                      */
    int pipe_cnt = 0;
    if (nbPipes) {
      /* Count number of '|' delimiters and detect the correct one */
      for (i = 0; seq.hdr[i] != '\0'; i++) {
        if (seq.hdr[i] == '|')
          pipe_cnt++;
        if (pipe_cnt == nbPipes) {
          /* Extract word after pipe delimiter  */
          char *s = &seq.hdr[i+1];
          char tmp[HDR_MAX];
          int j = 0;
          while (*s && !isspace(*s)) {
            tmp[j++] = *s++; 
          }
          tmp[j] = 0;
          /* Copy seq identifier to seq header  */
          strcpy (seq.hdr, tmp);
          break;
        }
      }
      /* If pipe_cnt = 0 leave header as it is  */
      if (pipe_cnt && pipe_cnt < nbPipes) { 
        strcpy(seq.hdr, "chrN");
      }  
    }
    if (options.debug)
      fprintf(stderr, "Sequence ID: %s\n", seq.hdr);
    /* Gobble sequence  */ 
    seq.len = 0;
    while ((res = fgets(buf, BUF_SIZE, input)) != NULL && buf[0] != '>') {
      char c;
      short int n;
      s = buf;
      while ((c = *s++) != 0) {
	if (isalpha(c)) {
	  c = (char) toupper(c);
	  switch (c) {
	  case 'A':
            n = 1;
            break;
	  case 'C':
            n = 2;
            break;
	  case 'G':
            n = 3;
            break;
	  case 'T':
            n = 4;
            break;
	  case 'N':
            n = 0;
	    break;
	  default:
            n = 0;
	    ;
	  }
          seq.len++;
	  if (seq.len >= mLen) {
	    mLen += BUF_SIZE;
	    seq.seq = realloc(seq.seq, (size_t)mLen * sizeof(short int));
	  }
	  seq.seq[seq.len] = n; /* sequence starts at seq.seq[1]  */
	}
      }
    }
    if (options.debug)
      fprintf(stderr, "Sequence length: %u\n", seq.len);
    /* We now have the (not nul terminated) sequence.
       Process it: on both or only forward directions   */
    if (seq.len != 0) {
      /* Scan the sequence for matches to the given PWM */
      if (options.forward) {
        if (wordLen == pwmLen) {
          scan_seq_1f(&seq);
        } else {
          scan_seq_2f(&seq);
        }
      } else { /* Scan both strands */
        if (wordLen == pwmLen) {
          scan_seq_1(&seq);
        } else {
          scan_seq_2(&seq);
        }
      }
    }
  }
  free(seq.hdr);
  free(seq.seq);
  if (input != stdin) {
    fclose(input);
  }
  return 0;
}


int
main(int argc, char *argv[])
{
  char *pwmFile = NULL;
  char *bgProb = NULL;
  char** tokens;
  int i = 0;

#ifdef DEBUG
  mcheck(NULL);
  mtrace();
#endif

  int option_index = 0;
  static struct option long_options[] =
      {
          /* These options may or may not set a flag.
             We distinguish them by their indices. */
          {"debug",   no_argument,       0, 'd'},
          {"help",    no_argument,       0, 'h'},
          {"coff",    required_argument, 0, 'c'},
          {"matrix",  required_argument, 0, 'm'},
          {"forward", no_argument,       0, 'f'},
          {"wordlen", required_argument, 0, 'i'},
          {"bgcomp",  required_argument, 0, 'b'},
          {"pipes",   required_argument, 0, 'n'},
          {"seqnorm", no_argument,       0, 'q'},
          {0, 0, 0, 0}
      };

  while (1) {
    int c = getopt_long(argc, argv, "dhfc:m:n:i:b:", long_options, &option_index); 
    if (c == -1)
      break;
    switch (c) {
    case 'd':
      options.debug = 1;
      break;
    case 'h':
      options.help = 1;
      break;
    case 'f':
      options.forward = 1;
      break;
    case 'c':
      cutOff = atoi(optarg);
      break;
    case 'm':
      pwmFile = optarg;
      break;
    case 'n':
      nbPipes = atoi(optarg);
      break;
    case 'i':
      wordLen = atoi(optarg);
      break;
    case 'b':
      bgProb = optarg;
      break;
    case '?':
      break;
    default:
      printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }
  if (optind > argc || pwmFile == NULL || cutOff == INT_MIN) {
    fprintf(stderr,
	    "Usage: %s [options] -m <pwm_file> -c <cut-off> [<] [< file_in] [> file_out]\n"
	    "      where options are:\n"
	    "        -d[--debug]            Print debug information\n"
	    "        -h[--help]             Show this help text\n"
	    "        -f[--forward]          Scan sequences in forward direction [def=bidirectional]\n"
	    "        -i[--wordlen] <len>    Length of the words in the word index array [def=%d]\n"
	    "        -b[--bgcomp]           Background model (residue priors), e.g. : 25,25,25,25\n"
	    "        -n[--pipes]            Number of pipe delimiters in FASTA header after which\n"
	    "                               The sequence identifier is expected to start [def=%d]\n"
	    "\n\tScan a DNA sequence file for matches to an INTEGER position weight matrix (PWM).\n"
            "\tThe DNA sequence file must be in FASTA format (<fasta_file>).\n"
            "\tThe matrix format is integer log-odds, where each column represents a nucleotide base\n"
            "\tin the following order: A, C, G, T. The program returns a list of matches in BED format.\n\n",
	    argv[0], wordLen, nbPipes);
    return 1;
  }
  /* Allocate space for both PWM and reverse PWM */
  pwm = (int **)calloc((size_t)pwmLen+1, sizeof(int *));   /* Allocate rows (PWM length + 1) */
  if (pwm == NULL) {
    fprintf(stderr, "Could not allocate matrix array: %s(%d)\n",
        strerror(errno), errno);
    return 1;
  }
  for (i = 0; i <= pwmLen; i++) {
    pwm[i] = calloc((size_t)NUCL, sizeof(int));                 /* Allocate columns (NUCL=5) */
    if (pwm[i] == NULL) {
      fprintf(stderr, "Out of memory\n");
      return 1;
    }
  }
  pwm_r = (int **)calloc((size_t)pwmLen+1, sizeof(int *)); /* Allocate rows (PWM length + 1) */
  if (pwm_r == NULL) {
    fprintf(stderr, "Could not allocate matrix array: %s(%d)\n",
        strerror(errno), errno);
    return 1;
  }
  for (i = 0; i <= pwmLen; i++) {
    pwm_r[i] = calloc((size_t)NUCL, sizeof(int));                  /* Allocate rows (NUCL=5) */
    if (pwm_r[i] == NULL) {
      fprintf(stderr, "Out of memory\n");
      return 1;
    }
  }
  /* Read Matrix from file */
  if ((pwmLen = read_pwm(pwmFile)) < 0)
    return 1;

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
  /* Word index length */
  if (wordLen > 15) {
    wordLen = 15;
  }
  if (wordLen <= 0) {
    wordLen = 7;
  }
  if (pwmLen < wordLen)
    wordLen = pwmLen;
  /* Read background model */
  if (bgProb != NULL) {
    tokens = str_split(bgProb, ',');
    if (tokens) {
      int i;
      for (i = 0; *(tokens + i); i++) {
        bgcomp[i] =  atof(*(tokens + i));
        free(*(tokens + i));
      }
      if (i != 4) {
        fprintf(stderr, "Number of TOKENS: %d\n", i);
        fprintf(stderr, "Please, set the bg base composition correctly: a comma-separated list of four numbers (e.g. 25,25,25,25)!\n");
        exit(1);
      }
      free(tokens);
    }
  }
  process_bgcomp(); 
  if (options.debug != 0) {
    if (fasta_in != stdin) {
      fprintf(stderr, "Fasta File : %s\n", argv[optind]);
    } else {
      fprintf(stderr, "Sequence File from STDIN\n");
    }
    fprintf(stderr, "Motif length: %d\n", pwmLen);
    fprintf(stderr, "Word index length: %d\n", wordLen);
    fprintf(stderr, "Weight Matrix: \n\n");
    for (int j = 1; j <= pwmLen; j++) {
      for ( int i = 1; i < NUCL; i++) {
        int mval = pwm[j][i];
        fprintf(stderr, " %7d ", mval);
      }
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
    if (!options.forward) {
    fprintf(stderr, "Reverse Weight Matrix:\n\n");
      for (int j = 1; j <= pwmLen; j++) {
        for ( int i = 1; i < NUCL; i++) {
          int mval = pwm_r[j][i];
          fprintf(stderr, " %7d ", mval);
        }
        fprintf(stderr, "\n");
      }
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "Background nucleotide frequencies:\n");
    for (int i = 1; i < NUCL; i++)
      fprintf(stderr, "bg[%i]=%1.2f ", i, bgcomp[i]);
    fprintf(stderr, "\n");
  }
  /* Re-scale matrices */
  process_pwm();
  if (pwmLen > wordLen)
    define_search_strategy();

  if (make_tables() != 0)
    return 1; 

  if (process_seq(fasta_in, argv[optind++]) != 0)
    return 1;

  /* Free PWMs structures */
  for (i = 0; i <= pwmLen; i++)
    free(pwm[i]);
  free(pwm);
  /* Free word index arrays (Scores and z link arrays) */
  free(z);
  free(ScoreF);
  free(Rfw);
  if (!options.forward) {
    free(ScoreR);
    free(Rrv);
    for (i = 0; i <= pwmLen; i++)
      free(pwm_r[i]);
    free(pwm_r);
  }
  return 0;
}
