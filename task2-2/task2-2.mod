param num_matrices 	>= 1, integer;	# Number of matrices in the data file to be read
param num_rows 		>= 1, integer;  # Number of rows
param num_cols 		>= 1, integer;  # Number of columns 

param min_dose >=0;    				   # Min dose for tumor
param max_dose >=0;    				   # Max dose for critical area
param min_var_bound >=0, <= min_dose;  # max variation for min_offset

param lambda >= 0, <= 1;

set BEAMS := 1 .. num_matrices;
set BROWS := 1 .. num_rows;
set BCOLS := 1 .. num_cols;


var S {b in BEAMS} >= 0; 			#strength of beam n

set MATS := 1 .. num_matrices; # set of matrices
set ROWS := 1 .. num_rows;	   # set of rows
set COLS := 1 .. num_cols;	   # set of columns

param matrix_value {MATS, ROWS, COLS} >= 0; # values for entries of each matrix
param tumor_value {ROWS,COLS} >= 0;
param crit_value {ROWS,COLS} >= 0;

set TUMOR := {i in ROWS, j in COLS: tumor_value[i,j] > 0}; 
set CRIT := {i in ROWS, j in COLS: crit_value[i,j] > 0}; 

var min_offset {(i,j) in TUMOR} >= 0, <= min_var_bound; 	# offset for min tumor dose
var max_offset {(i,j) in CRIT} >= 0;			  	# offset for max critical dose

# Pushing all variables to the maximum value of their corresponding indices
#minimize beamusage: sum {i in ROWS, j in COLS} sum{b in BEAMS}(S[b] * matrix_value[b,i,j]); #use this objective function to find minimum total dosage
/*maximize beamweight: lambda * sum {(i,j) in TUMOR}(sum{b in BEAMS}(S[b]*matrix_value[b,i,j])) -
	(1-lambda) * sum {(i,j) in CRIT}(sum{b in BEAMS}(S[b]*matrix_value[b,i,j])); #think this is right, but keep having infeasible solutions
*/
minimize Beam_and_Offsets: lambda * (sum {i in ROWS, j in COLS} sum{b in BEAMS}(S[b] * matrix_value[b,i,j]))
							+ (1 - lambda) * ((sum {(i,j) in TUMOR} min_offset[i,j]) + (sum {(i,j) in CRIT} max_offset[i,j]));

# Each variable at an index is >= to the maximum value at 
# the index across all matrices given.
subject to MinReq {(i,j) in TUMOR}: 
	sum{b in BEAMS} (S[b] * matrix_value[b,i,j]) >= min_dose - min_offset[i,j];

subject to MaxReq {(i,j) in CRIT}: 
	sum{b in BEAMS} (S[b] * matrix_value[b,i,j]) <= max_dose + max_offset[i,j];
