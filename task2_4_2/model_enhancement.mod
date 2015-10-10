param num_matrices >= 1, integer; # Number of matrices holding relative beam strengths
param num_beams >= 1, integer;    # Number of beams
param num_rows >= 1, integer;     # Number of rows
param num_cols >= 1, integer;     # Number of columns 
param num_lvl >=1, integer; 	  # Number of levels in the 3rd dimension

param min_dose >=0 ;    # Min dose for tumor
param max_dose >=0 ;     # Max dose for critical area

param lambda >=0, <=1 ;

set BEAMS := 0 .. num_beams-1; #beams;

var S{b in BEAMS} >=0; #strength of beam n

set MATS    := 1 .. num_matrices; # set of matrices
set ROWS    := 1 .. num_rows;	  # set of rows
set COLUMNS := 1 .. num_cols;	  # set of columns
set LVL     := 1 .. num_lvl;      # set of levels in the 3rd dimension

set TROWS := 1 .. num_rows;
set TCOLS := 1 .. num_cols;

set CROWS := 1 .. num_rows;
set CCOLS := 1 .. num_cols;


param matrix_value {MATS, ROWS, COLUMNS} >= 0; # values for entries of each matrix
param tumor_value {LVL,TROWS,TCOLS} >=0;
param crit_value {LVL,CROWS,CCOLS} >=0;

set TUMOR := {l in LVL, i in ROWS, j in COLUMNS: tumor_value[l,i,j]>0}; 
set CRIT := {l in LVL, i in ROWS, j in COLUMNS: crit_value[l,i,j]>0}; 

var max_offset {l in LVL, i in ROWS, j in COLUMNS} >= 0;
var min_offset {l in LVL, i in ROWS, j in COLUMNS} >=0, <=3;

# Pushing all variables to the maximum value of their corresponding indices


#minimize offs: sum{(l, i, j) in CRIT} max_offset[l, i, j]+ sum{(l,i,j) in TUMOR} min_offset[l,i,j];
minimize offs: sum{(l, i, j) in CRIT} max_offset[l, i, j];

# Each variable at an index is >= to the maximum value at 
# the index across all matrices given.

subject to MinReq {(l,i,j) in TUMOR}: 
	sum{b in BEAMS} (S[b] * matrix_value[(b*num_lvl+l),i,j]) >= min_dose;

subject to MaxReq {(l,i,j) in CRIT}: 
	sum{b in BEAMS} (S[b] * matrix_value[(b*num_lvl+l),i,j]) <= max_dose+max_offset[l,i,j];
