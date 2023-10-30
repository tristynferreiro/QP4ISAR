% A testing script used to measure the runtime of all RA and AF algorithm
% combinations considered for use in the QLP.
    % 
    % The CLI allows users to select which combination of RA and AF 
    % algorithms to test and saves the timing measurement for all 100 runs 
    % to a CSV file.

close all; clear all; clc; 
%% 1 Get input parameters
RA_selection = input("Range-alignment?\n 1 = Correlation\n 2 = Haywood\n");
AF_selection = input("Autofocus?\n 1 = Yuan \n 2 = Haywood\n");

%% (A) Setup path name for timing results
path_name ="";
if(RA_selection==1)
    path_name = "QLP_design_timing_CorrRA";
elseif(RA_selection==2)
    path_name = "QLP_design_timing_HaywoodRA";
end
if(AF_selection ==1)
    path_name = strcat(path_name,"_YuanAF");
elseif(AF_selection ==2)
    path_name = strcat(path_name,"_HaywoodAF");
end
path_name = strcat(path_name,".csv");
video_name = strcat(path_name,".avi");

%% (B) Time testing
time = zeros(100,1); % variable for storing time values
dataset = 'DAP_2010-10-14_09-43-33_010_zayaan_inbound_singlebin_P455_G1_sb_HRR.mat';
for i = 1:100
    fprintf("%d,",i); % Visual counter of timing progress
    t_start = tic;
    QLPDesignTimeTesting(128,50,RA_selection,AF_selection,video_name);
    time(i,1)=toc(t_start);
    close all;
end
%% (C) Save timing results to file
csvwrite(path_name, time); % write timing values to files
fprintf("\nDone.\n");