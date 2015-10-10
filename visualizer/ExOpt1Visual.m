%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Applied Math 121                          %
% Extreme Optimization 1: Let's cure cancer %
% Visualization script                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load beam output from AMPL
[file,path,success] = uigetfile({'*.out;*.txt'},...
                                'Please select the beam output file from AMPL, usually beam.out');
if success
    beam = importdata(strcat(path,file));
else
    error('EO1Visualizer:fileimport',...
          'Beam output file selection error. Please try again. If the problem persists, contact the course staff.');
end

% load tumor map
[file,path,success] = uigetfile({'*.out;*.txt'},...
                                'Please select the tumor map, usually tumor_raw.txt',...
                                path);
if success
    tumor = importdata(strcat(path,file));
else
    error('EO1Visualizer:fileimport',...
          'Tumor file selection error. Please try again. If the problem persists, contact the course staff.');
end

% load critical map
[file,path,success] = uigetfile({'*.out;*.txt'},...
                                'Please select the critical map, usually critical_raw.txt',...
                                path);
if success
    critical = importdata(strcat(path,file));
else
    error('EO1Visualizer:fileimport',...
          'Critical file selection error. Please try again. If the problem persists, contact the course staff.');
end

% check dimensions
if ~isequal(size(beam),size(tumor),size(critical))
    errormsg = strcat('Dimensions of beam, tumor, and critical region do not match. \nBeam:',num2str(size(beam,1)),'x',num2str(size(beam,2)),...
           '\nTumor:',num2str(size(tumor,1)),'x',num2str(size(tumor,2)),...
           '\nCritical:',num2str(size(critical,1)),'x',num2str(size(critical,2)));
    error('EO1Visualizer:sizemismatch',errormsg);
end
% check range of tumor and critical
if max(tumor(:))>1 || min(tumor(:))<0
    error('There are elements of tumor outside of the range [0,1]');
elseif max(critical(:))>1 || min(critical(:))<0
    error('There are elements of critical outside of the range [0,1]');
end
% round to ensure tumor and critical are binary
tumor = round(tumor);
critical = round(critical);

% get max value of beam output to normalize
beammax = max(beam(:));
% normalize beam
beamnorm = beam/beammax;

% compute complement of tumor as map of non-tumor areas
nontumor = imcomplement(tumor);

% create image
output = zeros(size(beam,1),size(beam,2),3);
output(:,:,1) = 0.5*tumor;
output(:,:,2) = beamnorm + 0.2*tumor;
output(:,:,3) = 0.2*nontumor + 0.3*beamnorm - 0.05*critical;
imshow(output);




