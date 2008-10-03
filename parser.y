%{
#include <stdio.h>
#include <stdlib.h>  
#include <string.h>
#include <search.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

#include "fm.h"

  extern int num_lines;
  extern polytope pgm;
%}
%union {
  double dvalue;
  long int livalue;
  int ivalue;
  char * svalue;
  term_node * tnp;
  constraint *constraint;
}
/* terminals */
%token <ivalue> EQ LE GE
%token <ivalue> PLUS MINUS
%token <svalue> VARIABLE
%token <dvalue> CONSTANT
%token <ivalue> INDEX

 /* nonterminals*/
%type <ivalue> sign
%type <ivalue> direction
%type <tnp> term 
%type <tnp> p1
%type <constraint> constraint file
%type <constraint> constraints
%%

file: constraints
;

constraints : constraint {
  insque(pgm.constraint_list=$1,NULL);
 }
| constraints constraint {
  insque($2,$1);
  $$=$2;
 }
;

direction : EQ | LE | GE
;

term : CONSTANT {$$=new_term_node("",$1,0);}
| VARIABLE{$$=new_term_node($1,1,0);}
| CONSTANT VARIABLE {$$=new_term_node($2,$1,0);}
;

p1 : sign term {$$=$2;$2->coefficient *= $1;} 
| sign term sign p1 {$$=$2;$2->coefficient *= $1; $4->coefficient*=$3;
   $2->next = $4;}
| term 
| term sign p1 {$3->coefficient *= $2 ; $1->next = $3;}
;

sign : PLUS | MINUS;

constraint : p1 direction CONSTANT {
  $$=new_constraint($1,$2,$3);
  /*printf("line: %d  ", num_lines);
    display_constraint($$);*/
  }
| p1 direction MINUS CONSTANT {
  $$=new_constraint($1,$2,-$4);
  /*printf("line: %d  ", num_lines);
    display_constraint($$);*/
  }
;
%%

