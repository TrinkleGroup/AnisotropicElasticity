#ifndef __DRAWFIG_H
#define __DRAWFIG_H

/*
  Program: drawfig.H
  Author:  D. Trinkle
  Date:    August 13, 2003
  Purpose: We define a class called drawfig that handles all of the
           output of a fig ver. 3.2 (man fig2dev for more info, or
	   check http://www.xfig.org/userman/fig-format.html ) file.

	   This class contains the following routines:

	   class drawfig
	     drawfig ()
	     drawfig (FILE)
	     drawfig (filename)
	       --constructors that open up the initial file only.
	         If no file is specified, stdout is assumed.

	     drawfig(FILE, portflag)
	     drawfig(FILE, portflag, a_scale, x_origin, y_origin)
	     drawfig(filename, portflag)
	     drawfig(filename, portflag, a_scale, x_origin, y_origin)
	       --construct and initialize.  portflag is true for
	         portrait mode, false for landscape.  If scale and
		 origin are present, they are used for converting
		 *double* numbers only.

	     init(portflag)
	       --see appropriate constructors

	     verbose()
	     verbose(VERBOSITY)
	       --returns, or sets, the level of verbosity.  I have
	         yet to decide what that does; I will probably use it
	         just to add comments of the untranslated coordinates
	         to the fig file.  In theory, we could go crazy and
	         add in text to specify the coordinates.  Initially
		 set to 0.

	     Each of the following returns, or sets, the appropriate
	     parameter:
	     
    	     pencolor() / pencolor(COLOR)
	     linethickness() / linethickness(THICKNESS)
	     linestyle() / linestyle(STYLE)
	     dotdist() / dotdist(DIST)
	     depth() / depth(DEPTH)

	     Then, we have output for lines:

	     line(x0, y0, x1, y1)
	       --if values are ints, then it is in xfig coord.  If they
	         are all doubles, it needs to be converted.

	     triangle(x0, y0, x1, y1, x2, y2)
	       --if values are ints, then it is in xfig coord.  If they
	         are all doubles, it needs to be converted.  Makes a filled
		 triangle using current fillstyles, etc.

	     vector(x, y, vx, vy)
	       --draw a vector (line with an arrow) starting at x, y
	         going in direction vx, vy

	     cvector(x, y, vx, vy)
	       --draw a centered vector (line with an arrow), centered
	         at x, y going in direction vx, vy

	     circle(x, y, r)
	       --draw a circle centered at x,y with radius r

	     text(x, y, s)
	       --puts string s at x,y

	     fillstyle(color, style) 
	       --sets the fill color and style (for circles).
		 
	     arrowstyle(type, filled, thick, length_percent, aspect)
	       --sets the type of arrows (for vectors).  The type is an
	         integer, filled is 0 / 1 for hollow vs. filled,
		 thickness is a double for the line thickness, and
		 the length is stored as a percentages of a vector's length,
		 while the aspect ratio means width = aspect*length.

	     textstyle(font, pointsize)
	       --sets the text font and pointsize.  We're a little sneaky,
	         though... if font = FONT_DEFAULT, we don't change it...
		 but we won't actually *output* to the fig the FONT_DEFAULT
		 tag.  Very sneaky, indeed.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <iomanip>
#include "io.H"

// ******************************* CONSTANTS ***************************

// Color definitions:
const int BLACK = 0;
const int WHITE = 7;
const int BLUE = 1;
const int GREEN = 2;
const int CYAN = 3;
const int RED = 4;
const int MAGENTA = 5;
const int YELLOW = 6;

// Fill options:
const int NOFILL = -1;
const int BLACKFILL = 0;
const int FULLFILL = 20;
const int WHITEFILL = 40;

// Font choices:
const int FONT_DEFAULT = -1;
const int FONT_TIMES_ROMAN = 0;
const int FONT_TIMES_ITALIC = 1;
const int FONT_TIMES_BOLD = 2;
const int FONT_TIMES_BOLD_ITALIC = 3;
const int FONT_AVANTGARDE_BOOK = 4;
const int FONT_AVANTGARDE_BOOK_OBLIQUE = 5;
const int FONT_AVANTGARDE_DEMI = 6;
const int FONT_AVANTGARDE_DEMI_OBLIQUE = 7;
const int FONT_BOOKMAN_LIGHT = 8;
const int FONT_BOOKMAN_LIGHT_ITALIC = 9;
const int FONT_BOOKMAN_DEMI = 10;
const int FONT_BOOKMAN_DEMI_ITALIC = 11;
const int FONT_COURIER = 12;
const int FONT_COURIER_OBLIQUE = 13;
const int FONT_COURIER_BOLD = 14;
const int FONT_COURIER_BOLD_OBLIQUE = 15;
const int FONT_HELVETICA = 16;
const int FONT_HELVETICA_OBLIQUE = 17;
const int FONT_HELVETICA_BOLD = 18;
const int FONT_HELVETICA_BOLD_OBLIQUE = 19;
const int FONT_HELVETICA_NARROW = 20;
const int FONT_HELVETICA_NARROW_OBLIQUE = 21;
const int FONT_HELVETICA_NARROW_BOLD = 22;
const int FONT_HELVETICA_NARROW_BOLD_OBLIQUE = 23;
const int FONT_NEW_CENTURY_SCHOOLBOOK_ROMAN = 24;
const int FONT_NEW_CENTURY_SCHOOLBOOK_ITALIC = 25;
const int FONT_NEW_CENTURY_SCHOOLBOOK_BOLD = 26;
const int FONT_NEW_CENTURY_SCHOOLBOOK_BOLD_ITALIC = 27;
const int FONT_PALATINO_ROMAN = 28;
const int FONT_PALATINO_ITALIC = 29;
const int FONT_PALATINO_BOLD = 30;
const int FONT_PALATINO_BOLD_ITALIC = 31;
const int FONT_SYMBOL = 32;
const int FONT_ZAPF_CHANCERY_MEDIUM_ITALIC = 33;
const int FONT_ZAPF_DINGBATS = 34;

const int FONT_MAX = 34;

// ******************************** CLASSES ****************************

class drawfig 
{
private:
  FILE* outfile;   // Output stream
  int FILEOPENED;  // flag for whether we opened the file or not
  int INITIALIZED; // flag for whether init's been called yet.  

  // Page definition:
  double a0;       // Scaling function
  double xc, yc;   // Origin shift
  int WIDTH, HEIGHT; // Width and height of page.
  int XC, YC;      // Center of page.
    
  // Parameters:
  int VERBOSE;
  int color_p;
  int thick_p;
  int style_p;
  double dotdist_p;
  int depth_p;

  // Arrow parameters:
  int arrow_type;
  int arrow_filled;
  double arrow_thick;
  double arrow_width;
  double arrow_height;

  // Fill parameters:
  int fillcolor_p;
  int fillstyle_p;

  // Text parameters:
  int font_p;
  double point_p;

  // This sets the default initial values for everything:
  void default_values () 
  {
    VERBOSE = 0;
    INITIALIZED = 0;
    color_p = BLACK;
    thick_p = 1;
    style_p = 0;
    dotdist_p = 0.;
    depth_p = 50;

    arrow_type = 2;
    arrow_filled = 1;
    arrow_thick = 1.;
    arrow_width = 0.1;
    arrow_height = 0.2;

    fillcolor_p = BLACK;
    fillstyle_p = NOFILL;
    
    font_p = FONT_COURIER;
    point_p = 10.;
  }

  // Get ready to print stuff out:
  void init(int portflag) 
  {
    if (INITIALIZED) {
      fprintf(stderr, "Tried to initialize drawfig twice?\n");    
    }
    else {
      fprintf(outfile, "#FIG 3.2\n");
      if (portflag) {
	fprintf(outfile, "Portrait\n");
	WIDTH = 10200;
	HEIGHT = 13200;
      }
      else {
	fprintf(outfile, "Landscape\n");
	WIDTH = 13200;
	HEIGHT = 10200;
      }
      fprintf(outfile, "Center\n");
      fprintf(outfile, "Inches\n");
      fprintf(outfile, "Letter\n");
      fprintf(outfile, "100.00\n");
      fprintf(outfile, "Single\n");
      fprintf(outfile, "-2\n");
      fprintf(outfile, "1200 2\n");
      INITIALIZED = -1;

      XC = WIDTH/2;
      YC = HEIGHT/2;
    }
  }

  // Convert from double to int for plotting.
  // Note: We want to make (0.0,0.0) at the center of the page,
  //   and positive y going *up* the page.
  inline int conv_x (double x) {return  (int)(a0*(x-xc)) + XC;}
  inline int conv_y (double y) {return -(int)(a0*(y-yc)) + YC;}

  // Convert *vectors*
  inline int conv_vx (double x) {return  (int)(a0*x);}
  inline int conv_vy (double y) {return -(int)(a0*y);}

public:
  // Constructors
  drawfig (FILE* file = stdout, int portflag = -1, 
	   double a_scale = 1.0, double x_origin = 0.0, double y_origin = 0.0)
    : outfile(file), a0(a_scale), xc(x_origin), yc(y_origin) {
    FILEOPENED = 0; default_values(); init(portflag);
  }

  drawfig (char* filename, int portflag = -1, 
	   double a_scale = 1.0, double x_origin = 0.0, double y_origin = 0.0)
    : a0(a_scale), xc(x_origin), yc(y_origin) {
    outfile = myopenw(filename); FILEOPENED = -1;
    if (outfile == NULL)
      fprintf(stderr, "Couldn't open file %s for output by drawfig.\n", filename);    
    default_values();
    init(portflag);
  }

  // Destructor
  ~drawfig() 
  {
    if (FILEOPENED)
      myclose(outfile);
  }

  // ********************* PARAMETER SETTING/READING *******************

  int verbose() const { return VERBOSE; }
  int verbose(int v) {VERBOSE = v; return VERBOSE;}

  int pencolor () const { return color_p; }
  int pencolor (int c) { color_p = c; return color_p; }
  
  int linethickness () const { return thick_p; }
  int linethickness (int t) { thick_p = t; return thick_p; }
  
  int linestyle () const { return style_p; }
  int linestyle (int s) { style_p = s; return style_p; }

  double dotdist () const { return dotdist_p; }
  double dotdist (double d) { dotdist_p = d; return dotdist_p; }
  
  int depth () const { return depth_p; }
  int depth (int d) { depth_p = d; return depth_p; }

  void arrowstyle(int a_type, int filled, double thick, 
		  double length_percent, double aspect) 
  {
    if ( (a_type >= 0) && (a_type <= 3) )
      arrow_type = a_type;
    if ( (filled == 0) || (filled == 1) )
      arrow_filled = filled;
    if ( thick > 0.)
      arrow_thick = thick;
    if ( length_percent > 0 ) {
      arrow_height = length_percent * 0.01;
      if ( aspect > 0 )
	arrow_width = length_percent * 0.01 * aspect;
    }
  }

  void fillstyle(int color, int style) 
  {
    if ( (color >= 0) && (color <= 31))
      fillcolor_p = color;
    if ( (style >= -1) && (style <= 40) )
      fillstyle_p = style;
  }

  void textstyle(int font, double pointsize) 
  {
    // if pointsize == 0, then no change.
    if (pointsize > 0.)
      point_p = pointsize;
    // Note: if font = FONT_DEFAULT == -1, we do nothing
    if ( (font_p >= 0) && (font_p <= FONT_MAX) )
      font_p = font;
  }

  // ************************** OUTPUT ROUTINES ************************

  // Make a line from (x0,y0) to (x1,y1)
  void line(int x0, int y0, int x1, int y1);
  void line(double x0, double y0, double x1, double y1);

  // Make a triangle from (x0,y0) to (x1,y1) to (x2,y2)
  void triangle(int x0, int y0, int x1, int y1, int x2, int y2);
  void triangle(double x0, double y0, double x1, double y1, double x2, double y2);

  // Make a vector at (x,y) in direction (vx,vy)
  void vector(int x, int y, int vx, int vy);
  void vector(double x, double y, double vx, double vy);

  // Make a centered vector at (x,y) in direction (vx,vy)
  void cvector(int x, int y, int vx, int vy);
  void cvector(double x, double y, double vx, double vy);

  // Make a circle centered at (x,y) with radius r
  void circle(int x, int y, int r);
  void circle(double x, double y, double r);

  void text(int x, int y, char* s);
  void text(double x, double y, char* s);
};


// *************************** OUTPUT ROUTINES *************************

void drawfig::line(int x0, int y0, int x1, int y1) 
{
  if ( ( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) ) < 1 ) return ;
  
  fprintf(outfile, "2 1 %d %d %d 0 %d 0 -1 %.3lf 0 0 0 0 0 2\n",
	  style_p, thick_p, color_p, depth_p, dotdist_p);
  fprintf(outfile, "%d %d %d %d\n", x0, y0, x1, y1);
}

void drawfig::line(double x0, double y0, double x1, double y1) 
{ line(conv_x(x0), conv_y(y0), conv_x(x1), conv_y(y1)); }


void drawfig::triangle(int x0, int y0, int x1, int y1, int x2, int y2) 
{
  if ( ( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) ) < 1 ) return ;
  
  fprintf(outfile, "2 3 %d %d %d %d %d 0 %d %.3lf 0 0 0 0 0 4\n",
	  style_p, thick_p, color_p, fillcolor_p, depth_p, fillstyle_p,
	  dotdist_p);
  fprintf(outfile, "%d %d %d %d %d %d %d %d\n", 
	  x0, y0, x1, y1, x2, y2, x0, y0);
}

void drawfig::triangle(double x0, double y0, double x1, double y1,
		       double x2, double y2) 
{ triangle(conv_x(x0), conv_y(y0), conv_x(x1), conv_y(y1), 
	   conv_x(x2), conv_y(y2)); }


void drawfig::vector(int x, int y, int vx, int vy)
{
  double len;
  len = sqrt(vx*vx + vy*vy);

  if ( (vx*vx + vy*vy) < 1 ) return;

  fprintf(outfile, "2 1 %d %d %d 0 %d 0 -1 %.3lf 0 0 0 1 0 2\n",
	  style_p, thick_p, color_p, depth_p, dotdist_p);
  fprintf(outfile, "%d %d %.3lf %.3lf %.3lf\n", arrow_type, arrow_filled,
	  arrow_thick, arrow_width*len, arrow_height*len);
  fprintf(outfile, "%d %d %d %d\n", x, y, x+vx, y+vy);
}

void drawfig::vector(double x, double y, double vx, double vy) 
{ vector(conv_x(x), conv_y(y), conv_vx(vx), conv_vy(vy)); }

 

void drawfig::cvector(int x, int y, int vx, int vy) 
{ vector(x - vx/2, y - vy/2, vx, vy); }

void drawfig::cvector(double x, double y, double vx, double vy) 
{ vector(x - vx*0.5, y - vy*0.5, vx, vy); }



void drawfig::circle(int x, int y, int r) 
{
  if ( r < 1 ) return;

  fprintf(outfile, "1 3 %d %d %d %d %d 0 %d %.3lf 1 0 ",
	  style_p, thick_p, color_p, fillcolor_p, depth_p, fillstyle_p,
	  dotdist_p);
  fprintf(outfile, "%d %d %d %d %d %d %d %d\n", x, y, r, r, 
	  x+r, y, x+r, y);
}

void drawfig::circle(double x, double y, double r) 
{ circle(conv_x(x), conv_y(y), conv_vx(r)); }



void drawfig::text(int x, int y, char* s) 
{
  double len, height;
  int i;

  for (i=0; s[i] != '\0'; ++i);
  len = 10.*point_p*i;
  height = 10.*point_p;

  fprintf(outfile, "4 1 %d %d 0 %d %.3lf 0. 6 %.3lf %.3lf %d %d %s\\001\n",
	  color_p, depth_p, font_p, point_p, height, len, x, y, s);
}


void drawfig::text(double x, double y, char* s) 
{ text(conv_x(x), conv_y(y), s); }

#endif
