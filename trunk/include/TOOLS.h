#ifndef TOOLSH
#define TOOLSH

extern void print_error(const char *str, BOOLEAN label);

extern void print_time(char printing, char *a);

extern BOOLEAN search_file(FILE *fp, const char *keyword);

extern BOOLEAN read_line(FILE * ffile, char * line, int MAXIMUM_LEN_FILE_LINE);

#endif

