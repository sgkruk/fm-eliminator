/* terminals */
%token <ivalue> EQ LE GE EOL ST
%token <ivalue> PLUS MINUS
%token <svalue> VARIABLE
%token <dvalue> CONSTANT
%token <ivalue> INDEX
%token <ivalue> MAXIMIZE MINIMIZE


%%
 /* Grammar for LP in an extension of CPLEX LP Format */

LP: Objective ST Constraints;


Objective: MaxOrMin LinearFunction;

MaxOrMin: MAXIMIZE | MINIMIZE;

LinearFunction: Term  SumOfTerms;


Term: CONSTANT | CONSTANT TIMES VARIABLE;

TIMES:  /* nothing, implied times */ | '*';

SumOfTerms: /* Nothing */ | Term SumOfTerms;



Constraints:  /* nothing */ | Constraint Constraints;

Constraint: LinearFunction RELATIONAL LinearFunction;

RELATIONAL: LE | EQ | GE;

%%
