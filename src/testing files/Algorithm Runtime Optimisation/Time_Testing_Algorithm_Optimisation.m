% A testing script used to measure the runtime of all RA and AF algorithms.
    % 
    % The CLI allows users to select which algorithm(s) to test and saves
    % the timing measurement for all 100 runs to a CSV file.

close all; clear all; clc; 

%% 1 Get input parameters
RA_selection = input("Range-alignment?\n 0 = none\n 1 = Correlation\n 2 = Haywood\n");
AF_selection = input("Autofocus?\n 0 = none\n 1 = Yuan \n 2 = Haywood\n");
optimised_selection = input("Optimised?\n 1 = Yes\n 2 = No\n");
%% 2 Load HRR profiles
HRR_profiles = csvread("Simulated_HRR_profiles.csv");
%% (A) Time Testing Setup
time = zeros(100,1); % variable for storing time values
%% 5 Range Alignment of Profiles 
if RA_selection ~=0
    % Range Align the HRR profiles using user selected algorithm
    ref_profile_number =1;
    if(RA_selection == 1) % CORRELATION RA
        for i = 1:100
            fprintf("%d,",i);% Visual counter of timing progress
            if(optimised_selection==1)
                path_name = strcat("Algorithm_Timing_CorrRA_optimised.csv");
                t_start = tic;
                RA_HRR_profiles = correlationRA(HRR_profiles,ref_profile_number);
                time(i,1)=toc(t_start);
            elseif(optimised_selection==2)
                path_name = strcat("Algorithm_Timing_CorrRA_unoptimised.csv");
                t_start = tic;
                RA_HRR_profiles = correlationRA_v0(HRR_profiles,ref_profile_number);
                time(i,1)=toc(t_start);
            end
            clear all; % Clean saved variables
        end
        csvwrite(path_name, time);% write timing values to files
        fprintf("\nDone.\n");
    elseif (RA_selection == 2) % HAYWOOD AF
        for i = 1:100
            fprintf("%d,",i);% Visual counter of timing progress
            if(optimised_selection==1)
                path_name = strcat("Algorithm_Timing_HaywoodRA_optimised.csv");
                t_start = tic;
                RA_HRR_profiles = HaywoodRA(HRR_profiles,ref_profile_number);
                time(i,1)=toc(t_start);
            elseif(optimised_selection==2)
                path_name = strcat("Algorithm_Timing_HaywoodRA_unoptimised.csv");
                t_start = tic;
                RA_HRR_profiles = HaywoodRA_v1(HRR_profiles,ref_profile_number);
                time(i,1)=toc(t_start);
            end
            clear all; % Clean saved variables
        end
        csvwrite(path_name, time);% write timing values to files
        fprintf("\nDone.\n");
    end
end

%% Autofocus of Profiles
if(AF_selection ~= 0)
    % Autofocus the HRR profiles using selected method
    if(RA_selection==0) % Use unaligned profiles if no RA selected.
        AF_RA_HRR_profiles = HRR_profiles; 
    end
    if(AF_selection == 1) % YUAN AF
        for i = 1:100
            fprintf("%d,",i);
            if(optimised_selection==1)
                path_name = strcat("Algorithm_Timing_YuanAF_optimised.csv");
                t_start = tic;
                AF_RA_HRR_profiles = YuanAF(RA_HRR_profiles);
                time(i,1)=toc(t_start);
            elseif(optimised_selection==2)
                path_name = strcat("Algorithm_Timing_YuanAF_unoptimised.csv");
                t_start = tic;
                AF_RA_HRR_profiles = YuanAF_v2(RA_HRR_profiles);
                time(i,1)=toc(t_start);
            end
            clear all; % Clean saved variables
        end
        csvwrite(path_name, time);% write timing values to files
        fprintf("\nDone.\n");
    elseif (AF_selection == 2) % HAYWOOD AF
        path_name = strcat("Algorithm_Timing_HaywoodAF_unoptimised.csv");
        for i = 1:100
            fprintf("%d,",i);
            if(optimised_selection==1)
                path_name = strcat("Algorithm_Timing_HaywoodAF_optimised.csv");
                t_start = tic;
                AF_RA_HRR_profiles = HaywoodAF(RA_HRR_profiles);
                time(i,1)=toc(t_start);
            elseif(optimised_selection==2)
                path_name = strcat("Algorithm_Timing_HaywoodAF_unoptimised.csv");
                t_start = tic;
                AF_RA_HRR_profiles = HaywoodAF_v1(RA_HRR_profiles);
                time(i,1)=toc(t_start);
            end
            clear all; % Clean saved variables
        end
        csvwrite(path_name, time);% write timing values to files
        fprintf("\nDone.\n");
    end
end