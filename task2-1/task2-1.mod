param num_matrices >= 1, integer; # Number of matrices in the data file to be read
param num_rows >= 1, integer;     # Number of rows
param num_cols >= 1, integer;     # Number of columns 

param min_dose >= 0;    # Min dose for tumor
param max_dose >= 0;    # Max dose for critical area

param lambda >= 0, <= 1;

set BEAMS := 1 .. num_matrices; #beams;
set BROWS := 1 .. num_rows;
set BCOLS := 1 .. num_cols;


var S {b in BEAMS} >= 0; #strength of beam n

set MATS    := 1 .. num_matrices; # set of matrices
set ROWS    := 1 .. num_rows;	  # set of rows
set COLUMNS := 1 .. num_cols;	  # set of columns


set TROWS := 1 .. num_rows;
set TCOLS := 1 .. num_cols;

set CROWS := 1 .. num_rows;
set CCOLS := 1 .. num_cols;


param matrix_value {MATS, ROWS, COLUMNS} >= 0; # values for entries of each matrix
param tumor_value {TROWS,TCOLS} >= 0;
param crit_value {CROWS,CCOLS} >= 0;

set TUMOR := {i in ROWS, j in COLUMNS: tumor_value[i,j] > 0}; 
set CRIT := {i in ROWS, j in COLUMNS: crit_value[i,j] > 0}; 

# Pushing all variables to the maximum value of their corresponding indices
minimize beamusage: sum {i in ROWS, j in COLUMNS} sum{b in BEAMS}(S[b] * matrix_value[b,i,j]); #use this objective function to find minimum total dosage
/* maximize beamweight: lambda * sum {(i,j) in TUMOR}(sum{b in BEAMS}(S[b]*matrix_value[b,i,j])) -
	(1-lambda) * sum {(i,j) in CRIT}(sum{b in BEAMS}(S[b]*matrix_value[b,i,j])); #think this is right, but keep having infeasible solutions
*/

# Each variable at an index is >= to the maximum value at 
# the index across all matrices given.
subject to MinReq {(i,j) in TUMOR}: 
	sum{b in BEAMS} (S[b] * matrix_value[b,i,j]) >= min_dose;

subject to MaxReq {(i,j) in CRIT}: 
	sum{b in BEAMS} (S[b] * matrix_value[b,i,j]) <= max_dose;
