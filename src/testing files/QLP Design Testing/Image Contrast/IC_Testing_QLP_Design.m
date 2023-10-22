% A testing script used to measure the image contrast (IC) of four frames 
% in the measured data set. 
%   All RA and AF algorithm combinations considered for use in the final 
%   QLP are tested. The values are saved to a CSV file.

close all; clear all; clc; 
%% (A) Setup path name for results
path_name ="QLP_design_IC_collated";
path_name = strcat(path_name,".csv");



%% (B) Image constrast testing
frame_center = [2464, 2970, 3827, 4189];
for frame = 1: size(frame_center,2)
   i=1; % keep track of number of columns
   for RA_selection = 1:2      %1 = Correlation, 2 = Haywood
       for AF_selection = 1:2  %1 = Yuan, 2 = Haywood
           [contrast_original(frame,1),contrast_after(frame,i)] = QLPDesignImageContrastTesting(128,frame_center(1,frame),RA_selection,AF_selection);
           i=i+1;
           fprintf(strcat("\nDone with ",num2str(RA_selection)," ",num2str(AF_selection),""));
           
       end
   end
end

%% (C) Save IC results to file
% Define the column names
column_names = {'Frame Center','Original','Focused CorrRA&YuanAF','Focused CorrRA&HaywoodAF','Focused HaywoodRA&YuanAF','Focused HaywoodRA&HaywoodAF'};
% Open the file for writing
file = fopen(path_name, 'w');
% Write the column names to the file
fprintf(file, '%s', column_names{1});
for i = 2:numel(column_names)
    fprintf(file, ',%s', column_names{i});
end
fprintf(file, '\n');
% Write the combined values to the file
combined_contrast = [frame_center',contrast_original,contrast_after];
dlmwrite(path_name, combined_contrast, '-append');
% Close the file
fclose(file);

% csvwrite(path_name, combined_contrast); % write timing values to files
fprintf("\nDone writing to file.\n");
