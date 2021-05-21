/*
#include <cstdlib>
#include "mycomplex.h"
#include "myconstants.h"
*/
#include <cstdio>
#include <ctime>
#include <cstring>
#include "mylogical.h"
#include "TOOLS.h"

using namespace std;

void print_time(char printing, char *a)

{

static clock_t start=0.0, end ;
static double total_time=0.0 ;

  double elapsed;

  end = clock();

  if ((double) end > (double) start) {
    elapsed = (double) (end - start) / CLOCKS_PER_SEC;
  } else {
    elapsed = (double) end / CLOCKS_PER_SEC;
  }

  total_time += elapsed;
  start = end;

  switch (printing) {
    case 0:
      printf("CPU time in %s: %f\n", a, elapsed);
    case 1:
      printf("Total CPU time: %f\n", total_time);
  }

  return;
}

/* this function reads entire line from an ASCII file */

BOOLEAN read_line(FILE * ffile, char * line, int MAXIMUM_LEN_FILE_LINE) {
  char c;
  int len;
  BOOLEAN line_too_long;

  len = 0;

  line[0] = '\0';

  if (feof(ffile) != 0) {
    return FALSE;
  }

  line_too_long = FALSE;

  while (fscanf(ffile, "%c", &c) == 1) {

    if (c == 0x0A) { /* !!! UNIX only !!! */

      break;
    }
    if (len >= MAXIMUM_LEN_FILE_LINE) {

      line_too_long = TRUE;
    } else {

      line[len] = c;
      len++;
      line[len] = '\0';
    }
  }

  if (line_too_long == TRUE) {

    fprintf(stderr, "ERROR: The line is too long. The limit is %d characters.\n", MAXIMUM_LEN_FILE_LINE);
  }

  return TRUE;
}

BOOLEAN search_file(FILE *fp, const char *keyword)

{

  BOOLEAN err;
  char ch[6];
  int i;

  /******************************************************
   * search_file START                                  *
   * CARE: word may not be longer than 6 characters     *
   ******************************************************/

  err = FALSE;
  while ((fscanf(fp, "%6s", ch)) != EOF) {
    if (!strcmp(ch, keyword)) {
      //fprintf(file_out, " Keyword, %6s found \n", ch );
      err = TRUE;
      break;
    }
  }

  return err;
}

BOOLEAN search_file_long_word(FILE *fp, char *keyword)

{

  BOOLEAN err;
  char ch[12];

  /*******************************************************
   * search_long_word START                              *
   * CARE: word may not be longer than 12 characters *
   ******************************************************/

  err = FALSE;
  while ((fscanf(fp, "%12s", ch)) != EOF) {
    if (!strcmp(ch, keyword)) {
      //fprintf(file_out, " Keyword, %6s found \n", ch );
      err = TRUE;
      break;
    }
  }

  return err;
}

