struct term_node {
  struct term_node *next;
  struct term_node *previous;
  char *variable;
  double coefficient;
  unsigned short column_index;
  int ordinal;
};
struct constraint {
  struct constraint *next;
  struct constraint *previous;
  struct term_node *first;
  short direction; /* -1 <=    0 ==      +1 >= */
  double rhs;
};

struct variable_node{
  struct variable_node *next;
  struct variable_node *previous;
  char *varname;
  unsigned short column_index;		
  unsigned char keep;		/* or eliminate */
};

typedef struct variable_node variable_node;
typedef struct constraint constraint;
typedef struct term_node term_node;

struct row_node {
  struct row_node *next;
  struct row_node *previous;
  unsigned short deleterow;
  short direction;
  gsl_vector *v;
};
typedef struct row_node row_node;

struct polytope {
  constraint *constraint_list;
  variable_node *variable_list;
  variable_node **variable_array;
  unsigned int nb_constraints;
  unsigned int nb_variables;
  row_node *row_list;
  unsigned int nb_rows;
  unsigned int nb_to_eliminate;
};
typedef struct polytope polytope;


row_node *new_row_node(int size);
term_node *new_term_node(char *variable,double coefficient,int column_index);
void display_term_node(term_node *tn);
void display_linear_function(term_node *tn);
constraint *new_constraint(term_node *tn,short direction,double rhs);
void display_constraint(constraint *c);
void extract_vars(polytope *pgm);
void display_variable_node(variable_node *vn);
variable_node *new_variable_node(char *varname);
