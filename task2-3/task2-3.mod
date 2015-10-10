param num_matrices 	>= 1, integer; 	# Number of matrices in the data file to be read
param num_rows 		>= 1, integer;  # Number of rows
param num_cols 		>= 1, integer;  # Number of columns 

param min_dose >=0 ;	# Min dose for tumor
param max_dose >=0 ;	# Max dose for critical area

param lambda_t >=0, <=1;
param lambda_c >=0, <=1;
param lambda_b >=0, <=1;
param lambda_o >=0, <=1;
param width >= 0;
set bound := -width .. width;

set BEAMS := 1 .. num_matrices;
set BROWS := 1 .. num_rows;
set BCOLS := 1 .. num_cols;

var S{b in BEAMS} >=0; 				# strength of beam n
var ma >= 0, <= min_dose; 	# offset for min tumor dose
var mb >= 0;			  	# offset for max critical dose
var maa >=0;
var mab >= 0;

set MATS := 1 .. num_matrices; 	# set of matrices
set ROWS := 1 .. num_rows;	  	# set of rows
set COLS := 1 .. num_cols;	  	# set of columns

param matrix_value {MATS, ROWS, COLS} >= 0; # values for entries of each matrix
param tumor_value {ROWS,COLS} >=0;
param crit_value {ROWS,COLS} >=0;

set TUMOR  := {i in ROWS, j in COLS: tumor_value[i,j]>0}; 
set CRIT   := {i in ROWS, j in COLS: crit_value[i,j]>0}; 
set BORDER_CALC := {i in ROWS, j in COLS: exists {k in bound, l in bound}
				(i + k >= 1) and (i + k <= num_rows) and (j + l >= 1) and (j + l <= num_cols)
				and (i + k, j + l) in CRIT};
set BORDER := BORDER_CALC diff CRIT diff TUMOR;
set OTHER := {i in ROWS, j in COLS} diff TUMOR diff CRIT diff BORDER;

# Pushing all variables to the maximum value of their corresponding indices
/*
minimize Beam_and_Offsets: lambda_t * (sum {(i,j) in TUMOR} sum{b in BEAMS}(S[b] * matrix_value[b,i,j]))
							+ lambda_c * (sum {(i,j) in CRIT} sum{b in BEAMS}(S[b] * matrix_value[b,i,j]))
							+ lambda_b * (sum {(i,j) in BORDER} sum{b in BEAMS}(S[b] * matrix_value[b,i,j]))
							+ lambda_o * (sum {(i,j) in OTHER} sum{b in BEAMS}(S[b] * matrix_value[b,i,j]));
*/
/*
maximize radiation:
							+ lambda_b * (sum {(i,j) in BORDER} sum{b in BEAMS}(S[b] * matrix_value[b,i,j]))
		
							+ (1 - lambda_t - lambda_c - lambda_b - lambda_o) * (min_offset + max_offset);
*/

maximize fun: sum{(i,j) in TUMOR} (sum{b in BEAMS} (S[b] * matrix_value[b,i,j]) -10); 

# Each variable at an index is >= to the maximum value at 
# the index across all matrices given.

/*
subject to MinReq {(i,j) in TUMOR}: 
	sum{b in BEAMS} (S[b] * matrix_value[b,i,j]) = min_dose - ma + mb;
*/

subject to MaxReq {(i,j) in CRIT}: 
	sum{b in BEAMS} (S[b] * matrix_value[b,i,j]) <= max_dose - mab;
