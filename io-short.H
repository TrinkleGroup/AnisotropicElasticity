#ifndef __IO_H
#define __IO_H

/*
  Program: io-short.H
  Author:  D. Trinkle
  Date:    2008 May 7
  Purpose: Short and sweet io routines
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <libgen.h>

// *************************** I/O SUBROUTINES ***************************
inline FILE* myopenr(char *filename) 
{ if (filename[0] == '-') return stdin;
  else                    return fopen(filename, "r"); }

inline FILE* myopenw(char *filename) 
{ if (filename[0] == '-') return stdout;
  else                    return fopen(filename, "w"); }
  
inline void myclose(FILE *f) 
{ if ((f != stdin) && (f != stdout) && (f != NULL)) fclose(f); }

inline int has_error (int ERROR, int ERROR_CODE) 
{ return ( (ERROR & ERROR_CODE) != 0 ); }

inline int is_blank (char c) {return (c == ' ') || (c == '\t') || (c == '\n');}

inline void nextnonblank (char* &p) 
{ for ( ; ((*p)!='\0') && is_blank(*p); ++p) ;}

inline void nextblank (char* &p) 
{ for ( ; ((*p)!='\0') && (!is_blank(*p)); ++p) ;}

const char COMMENT_CHAR = '#';
const char EOF_CHAR = '&';

inline void nextnoncomment (char* dump, const int &size, FILE* infile) 
{ do {char *p_fgets=fgets(dump, size, infile);} 
  while ((!feof(infile)) && (dump[0] == COMMENT_CHAR)); }

#endif
