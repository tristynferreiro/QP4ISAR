%close all; clear all; clc; 

%% Load CSIR dataset 

load('DAP_2010-10-14_09-43-33_010_zayaan_inbound_singlebin_P455_G1_sb_HRR.mat');
HRRProfilesAll = sb_HRR.G1.HRR_NoMC_calib.';
ProfileRepetitionFreq =  1/sb_HRR.G1.pattern_time; 
NumRangeBins = size(HRRProfilesAll,1);
NumOfProfiles = size(HRRProfilesAll,2);
Range_axis = sb_HRR.G1.xaxis_downrange_m; 

% Select a subset of profiles
MiddleProfile = 2464;                                              % 2464, 4189, 2970, 3827
CPTWL = 128;   
StartProfile = MiddleProfile - CPTWL/2;                           
StopProfile = MiddleProfile + CPTWL/2 - 1;
ProfilesToProcess = StopProfile - StartProfile;
DopplerAxis_Hz = (-ProfilesToProcess/2:1:ProfilesToProcess/2-1)*ProfileRepetitionFreq/ProfilesToProcess;

HRR_profiles = circshift(HRRProfilesAll(StartProfile:StopProfile, :), [0 50 ]);

test_prof =HRRProfilesAll(StartProfile:StopProfile, :);
% figure; imagesc(Range_axis, 1:size(HRR_profiles,1), 20*log10(abs(HRR_profiles)));
% figure; imagesc(Range_axis, 1:size(test_prof,1), 20*log10(abs(test_prof)));

%% Range Alignment of Profiles 
% Range Align the HRR profiles using correlation method
ref_profile_number =1;

% f = @() correlationRA(HRR_profiles,ref_profile_number); % handle to function
% timeit(f)
% f = @() HaywoodRA(HRR_profiles, ref_profile_number); % handle to function
% timeit(f)

%[RA_HRR_profiles_hay] = HaywoodRA(HRR_profiles,ref_profile_number);
tic
[RA_HRR_profiles_corr] = HaywoodRA(HRR_profiles,ref_profile_number);
%figure; imagesc(Range_axis, 1:size(RA_HRR_profiles_corr,1), 20*log10(abs(RA_HRR_profiles_corr)));
toc
%RA_HRR_profiles_haywood = HaywoodRA(HRR_profiles, ref_profile_number);
%% Autofocus of Profiles
% Apply AF to the RA HRR profiles using

% f = @() HaywoodAF(RA_HRR_profiles_corr); % handle to function
% timeit(f)
tic
[AF_RA_HRR_profiles_corr] = HaywoodAF(RA_HRR_profiles_corr);
toc

% f = @() HaywoodAF(RA_HRR_profiles_corr); % handle to function
% timeit(f)



