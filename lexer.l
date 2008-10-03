%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <search.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include "fm.h"
#include "parser.h"
extern int num_lines;
extern int num_chars;
%}
char [A-Za-z]
num [0-9,]
eq "="
ge ">="
le "<="
plus "+"
minus "-"
variable {char}+{num}*{arrayindex}?
constant {num}+
arrayindex "("{num}+")"
%%
\n ++num_lines;
{variable} { yylval.svalue = strdup(yytext);
             return VARIABLE; }
{eq} { yylval.ivalue=0; return EQ; }
{le} { yylval.ivalue=-1; return LE; }
{ge} { yylval.ivalue=1; return GE; }
{plus} { yylval.ivalue=1; return PLUS;}
{minus} { yylval.ivalue=(-1); return MINUS;}

{constant} { yylval.dvalue = atof(yytext);
            return CONSTANT; }
%%
int yywrap()
{
return 1;
}