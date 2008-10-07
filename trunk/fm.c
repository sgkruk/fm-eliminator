#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <search.h>
#include <unistd.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

#include "setoper.h"
#include "cdd.h"
#include "math.h"

#include "fm.h"

int num_lines=0;
int num_chars=0;
polytope pgm;
extern FILE *yyin;
extern int yyparse(void);

/*
  The code related to the reading of the constraints and storing
  in a general-purpose structure.
 */
term_node *new_term_node(char *variable,double coefficient,int column_index){
  term_node *tnp;
  tnp = (term_node *)malloc(sizeof(term_node));
  tnp->variable=strdup(variable);
  tnp->column_index=column_index;
  tnp->coefficient=coefficient;
  tnp->next=0;

  return tnp;
}
void display_term_node(term_node *tn){
  printf(" %+f %s [%d]", tn->coefficient, tn->variable, tn->column_index);
}
void display_linear_function(term_node *tn){
  term_node *c;
  for (c = tn; c; c = c->next){
    display_term_node(c);
  }
}
void display_polytope(polytope p){
  constraint *c=0;
  variable_node *vn=0;
  int i;
  printf("Polytope is determined by %d constraints and %d variables\n",
	 p.nb_constraints,p.nb_variables);
  for (c = p.constraint_list; c; c = c->next){
    display_constraint(c);
  }
  for(i=0,vn = p.variable_list; vn; vn=vn->next,i++){
    if (vn!=p.variable_array[i])
      printf("  Structure inconsistent !!!");
    /*display_variable_node(vn);*/
  }
}
void display_variable_node(variable_node *vn){
  printf("%s at column %d   (keep %d)\n", 
	 vn->varname, 
	 vn->column_index,
	 vn->keep);
}
constraint *new_constraint(term_node *tn,short direction,double rhs){
  constraint *c;
  c=(constraint *)malloc(sizeof(constraint));
  c->next=0;
  c->direction=direction;
  c->rhs=rhs;
  c->first=tn;

  return c;
}
void display_constraint(constraint *c){
  display_linear_function(c->first);
  switch(c->direction){
  case -1: printf("<=");break;
  case 0: printf("==");break;
  case 1: printf(">=");break;
  }
  printf(" %f\n",c->rhs);
}
variable_node *valid_variable(char *varname){
  ENTRY hashelement;
  ENTRY *found;
  hashelement.key=varname;
  if ((found=hsearch(hashelement,FIND)))
    return (variable_node *)found->data;
  else
    return 0;
}
variable_node *valid_regex_variable(polytope p,char *varname){
  variable_node *vn=0;
  for(vn = p.variable_list; vn; vn=vn->next){
	if (vn->keep) 
	  continue;					/* Been there already */
	else if (!strncmp(vn->varname,varname,strlen(varname)))
	  return vn;
  }
  return 0;
}
void extract_vars(polytope *p){
  constraint *c;
  term_node *tn;
  variable_node *vn,*lvn;
  unsigned short column_index=0;	/* Column index of var */
  ENTRY hashelement;
  ENTRY *found;
  unsigned int nb_constraints=0;
  unsigned int nb_variables=0;
  int i;
  /* setup the hash table */
  if (!hcreate(1000)){		/* Asking for trouble here */
    printf("No more memory, forget it\n");
    exit(-1);
  }
  p->nb_to_eliminate=0;
  for (c = p->constraint_list; c; c = c->next){
    nb_constraints++;
    for (tn=c->first; tn; tn=tn->next){
      hashelement.key=strdup(tn->variable);
      if (!(found=hsearch(hashelement,FIND))){		/* Not there yet */
		nb_variables++;
		vn=new_variable_node(tn->variable);
		if (!vn->keep) p->nb_to_eliminate++;
		if (!p->variable_list){
		  p->variable_list=vn;
		  insque(vn,0);
		} else {
		  insque(vn,lvn);
		}
		vn->column_index=column_index++;		/* set column index */
		tn->column_index=vn->column_index;
		lvn=vn;
		hashelement.data=vn;
		hsearch(hashelement,ENTER);
      } else {			/* Already in there */
		vn=(variable_node *)found->data;
		tn->column_index=vn->column_index;
      }
    }
  }
  p->nb_constraints=nb_constraints;
  p->nb_variables=nb_variables;
  /* Now the reverse index on the variable names */
  p->variable_array=
	(variable_node **)malloc(nb_variables*sizeof(variable_node *));
  for(i=0,vn=p->variable_list; vn; vn=vn->next,i++)
    p->variable_array[i]=vn;
}
/*
  The following code handles the numerical part of the Fourier-Motzkin 
  elimination process.
 */
row_node *new_row_node(int size){
  row_node *rn = (row_node *)malloc(sizeof(row_node));
  rn->v = gsl_vector_alloc(size);

  return rn;
}
row_node *new_vector_from_constraint(constraint *c,int vsize){
  term_node *tn;
  row_node *rn;
  rn=new_row_node(vsize);
  rn->direction=c->direction;
  for(tn=c->first; tn; tn=tn->next){
    gsl_vector_set(rn->v,tn->column_index,(tn->coefficient));
  }
  gsl_vector_set(rn->v,vsize-1,(c->rhs));
  /*display_constraint(c);
	gsl_vector_fprintf(stdout,rn->v,"%f");*/

  return rn;
}

void fill_rows_from_constraints(polytope *p){
  int nb_variables = p->nb_variables;
  constraint *c;
  row_node *rn,*rn0,*lrn;
  p->nb_rows=0;
  for (c = p->constraint_list; c; c = c->next){
	switch (c->direction){
	case 0: 
	case 1: 
	case -1:
	  rn = new_vector_from_constraint(c,nb_variables+1);
	  rn0=0; p->nb_rows++; break;
	case 2:
	  rn = new_vector_from_constraint(c,nb_variables+1);
	  rn0 = new_vector_from_constraint(c,nb_variables+1);
	  p->nb_rows+=2; 
	  break;
	}
    if (!p->row_list){
      p->row_list=rn;
      insque(rn,0);
    } else {
      insque(rn,lrn);
    }
    lrn=rn;
	if (rn0){
	  insque(rn0,rn);
	  lrn=rn0;
	}
  }
}
variable_node *new_variable_node(char *varname){
  variable_node *vn = (variable_node *)malloc(sizeof(variable_node));
  vn->varname=strdup(varname);
  return vn;
}
void display_rows(polytope p){
  row_node *rn;
  unsigned int i;
  variable_node *vn;
  for(rn=p.row_list; rn; rn=rn->next){
    for(i=0;i<p.nb_variables;i++){
      double val=gsl_vector_get(rn->v,i);
      if (val){
		vn=p.variable_array[i];
		printf("%+5.5f %4s ",val,vn->varname);
      } else {
		printf("              ");
      }
    }
	switch(rn->direction){
	case -1: printf(" <= "); break;
	case  0: printf("  = "); break;
	case +1: printf(" >= "); break;
	}
    printf(" %5.f\n", gsl_vector_get(rn->v,p.nb_variables));
	
  }
}
/* Add the new row in the structure (at end OUCH!) */
void fm_eliminate_pair(polytope *p,row_node *rp,row_node *rn){
  row_node *r = new_row_node(p->nb_variables+1);
  row_node *last;
  gsl_vector_memcpy(r->v,rn->v);
  gsl_vector_add(r->v,rp->v);
  for(last=p->row_list;last && last->next;last=last->next)
	/* nothing*/;
  insque(r,last);
  return;
}
/* Delete all rows that were marked delete and adjust structure p */
int delete_processed_rows(polytope *p){
  row_node *r;
  for (r=p->row_list;r;r=r->next)
	if (r->deleterow){
	  if (r==p->row_list)
		p->row_list=r->next;	/* In case we are deleting the first */
	  remque(r);
	}
  
  return 0;
}
/* Now we do the elimination */
void fm(polytope *p){
  int verbose=0;
  row_node *rp,*rn;
  int iv;
  variable_node *vn;
  int nb_constraints = p->nb_constraints;
  int nb_variables = p->nb_variables;
  for (iv=0; iv<nb_variables; iv++){
    /* for each variable, check if we keep it or not and act */
    vn=p->variable_array[iv];
    if ((vn->keep)) 
      continue;
    else {
      int 	  i=vn->column_index;
      if (verbose) printf("We will eliminate %s at position %d\n", 
			  vn->varname,
			  vn->column_index);
      for(rp=p->row_list;rp;rp=rp->next){
	double pv;							 
	/* Positive  */
	if ((pv=gsl_vector_get(rp->v,i))>0){ 
	  /* Found a positive coefficient */
	  rp->deleterow=1;						 
	  /* To be deleted at end */
	  if (pv != 1.0)
	    gsl_vector_scale(rp->v,1/pv);
	  if (verbose) display_rows(*p);
	  for(rn=p->row_list;rn;rn=rn->next){
	    double nv;							 
	    /* Negative */
	    if ((nv=gsl_vector_get(rn->v,i))<0){ 
	      /* Negative coefficient */
	      rn->deleterow=1;						 
	      /* To be deleted at end */
	      if (nv != -1.0)
		gsl_vector_scale(rn->v,-1/nv);
	      if (verbose) printf("positive %f negative %f\n",pv,nv);
	      /* Do the work */
	      fm_eliminate_pair(p,rp,rn);
	    }
	  }		  
	}
      }
      nb_constraints=delete_processed_rows(p);
    }
  }
}
/* Try to use CDDLIB to do that */
void project(polytope *p,FILE *outfile,int formatcdd){
  row_node *rn;
  variable_node *vn;
  double val;

  dd_MatrixPtr M=NULL,M1;  

  dd_colset delset;
  dd_colrange d;
  dd_rowindex newpos;
  dd_ErrorType err=dd_NoError;
  dd_rowset redset,impl_linset;

  unsigned i,j;

  dd_set_global_constants();  /* First, this must be called. */

  M=dd_CreateMatrix(p->nb_rows,p->nb_variables+1);
  M->representation=dd_Inequality;
  M->numbtype=dd_GetNumberType("real");
  
  for(i=0,rn=p->row_list; rn; i++,rn=rn->next){
    double r;
    if (rn->direction<0)
      r=-1.0;
    else if (rn->direction>0)
      r=1.0;
    else {
      r=-1.0;
      set_addelem(M->linset,i+1);
    }
    for(j=0;j<p->nb_variables;j++){
      dd_set_d(M->matrix[i][j+1],(r)*gsl_vector_get(rn->v,j));
    }
    dd_set_d(M->matrix[i][0],gsl_vector_get(rn->v,p->nb_variables));
  }
  /*dd_WriteMatrix(outfile, M);*/
  /* need to set the variables to eliminate */
  d=p->nb_variables+1;
  set_initialize(&delset,d);
  for (i=0; i<p->nb_variables; i++){
    /* for each variable, check if we keep it or not and act */
    vn=p->variable_array[i];
    if (!(vn->keep)) {
      long t=i+2;
      set_addelem(delset,t);
    }
  }
  M1=dd_BlockElimination(M, delset, &err);
  if (err!=dd_NoError) dd_WriteErrorMessages(stderr,err);
  dd_MatrixCanonicalize(&M1,&impl_linset,&redset,&newpos,&err);
  
  if (formatcdd)
    dd_WriteMatrix(outfile, M1);
  else {
    /* Print it for now */
    for(i=0; i<M1->rowsize; i++){
      vn=p->variable_list;
      for(j=1; j<M1->colsize; j++){
	val=(-1.0)*dd_get_d(M1->matrix[i][j]);
	while(!vn->keep)
	  vn=vn->next;
	if (val){
	  fprintf(outfile,"%+5.5f ",val );
	  fprintf(outfile,"%-4s ", vn->varname);
	} else {
	  fprintf(outfile,"              ");
	}
	vn=vn->next;
      }
      val=dd_get_d(M1->matrix[i][0]);
      if (set_member(i+1,M1->linset))
	fprintf(outfile,"== %5.5f \n",val );
      else
	fprintf(outfile,"<= %5.5f \n",val );
    }
  } /* end of print */
  dd_FreeMatrix(M);
  dd_FreeMatrix(M1);
}			  

/*
  Mainline 
*/
void usage(char *pgmname){
  printf("Usage is %s -i infile -o outfile (-[pr] var)+ [-f [0|1]]\n",pgmname);
  printf("Order of the parameters is important!\n");
  printf("-p var (exact name) or -r var (beginning of name)\n");
  printf("-f 1 for cdd format (-f 0 is default LP format)\n");
  
  return;
}
int main(int argc, char **argv)
{
  FILE *outfile=stdout;		/* Default is stdout */
  int verbose=0;
  int switchval;
  int formatcdd=0;		/* Format: Fukuda=1 LP=0 */
  variable_node *vn;
  while((switchval=getopt(argc, argv, "hvr:p:i:o:f:")) != -1){
    switch(switchval){
    case 'h':
      usage(argv[0]);
      exit(0);
    case 'i':
      if ((yyin=fopen(optarg,"r"))==0){
	perror(optarg);
	exit(0);
      }
      yyparse();
      extract_vars(&pgm);
      break;
    case 'o':
      if ((outfile=fopen(optarg,"w"))==0){
	perror(optarg);
	exit(0);
      }
      break;
    case 'f':
      formatcdd=atoi(optarg);
      break;
    case 'v': verbose++;
      display_polytope(pgm);  
      break;
    case 'p': 
      if ((vn=valid_variable(optarg))){
	if (verbose) printf("Projecting on %s\n",vn->varname);
	vn->keep=1;
      } else {
	if (verbose) printf("variable %s is not in use!\n",optarg);
      }
    case 'r':
      while((vn=valid_regex_variable(pgm,optarg))){
	if (verbose) printf("Projecting on %s\n",vn->varname);
	vn->keep=1;
      }
      break;
    default:
      usage(argv[0]);
      exit(0);
    }
  }
  
  fill_rows_from_constraints(&pgm);
  if (verbose>3) {
	printf("The rows before elimination are\n");
	display_rows(pgm);
  }
  /*
  fm(&pgm);
  printf("The rows after elimination are\n");
  display_rows(pgm);
  return 0;
  */
  project(&pgm,outfile,formatcdd);				/* Using CDD */
  return 0;
}
int yyerror(char *msg)
{
  printf("Error at line %d : %s \n", num_lines,msg);
  return 0;
}
