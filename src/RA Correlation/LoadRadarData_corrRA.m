close all; 
clear all;
clc; 

%% Load CSIR dataset 

load('DAP_2010-10-14_09-43-33_010_zayaan_inbound_singlebin_P455_G1_sb_HRR.mat');
HRRProfilesAll = sb_HRR.G1.HRR_NoMC_calib.';
ProfileRepetitionFreq =  1/sb_HRR.G1.pattern_time; 
NumRangeBins = size(HRRProfilesAll,1);
NumOfProfiles = size(HRRProfilesAll,2);
Range_axis = sb_HRR.G1.xaxis_downrange_m; 

% Plot the HRR profiles
figure; imagesc(sb_HRR.G1.xaxis_downrange_m, 1:size(HRRProfilesAll,1), 20*log10(abs(HRRProfilesAll))); 
xlabel('Range (m)');
ylabel('Profile Number');
title('All HRR profiles');
colormap('jet');
colorbar;

% Select a subset of profiles
MiddleProfile = 2464;                                              % 2464, 4189, 2970, 3827
CPTWL = 128;   
StartProfile = MiddleProfile - CPTWL/2;                           
StopProfile = MiddleProfile + CPTWL/2 - 1;
ProfilesToProcess = StopProfile - StartProfile;
DopplerAxis_Hz = (-ProfilesToProcess/2:1:ProfilesToProcess/2-1)*ProfileRepetitionFreq/ProfilesToProcess;

% Plot Unaligned HRR Profiles
HRR_profiles = circshift(HRRProfilesAll(StartProfile:StopProfile, :), [0 50 ]);
figure; imagesc(sb_HRR.G1.xaxis_downrange_m, 1:size(HRR_profiles,1), 20*log10(abs(HRR_profiles))); 
xlabel('Range (m)');
ylabel('Profile Number');
title('Subset of HRR profiles');
colormap('jet');
colorbar;

%% Plot unfocused ISAR image - no range alignment, no autofocus
% apply hamming window function
window = hamming(size(ProfilesToProcess,1));

ISAR_image = fftshift(fft(HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension                    

ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 

% Plot unfocused ISAR image
figure; imagesc( Range_axis, DopplerAxis_Hz, ISAR_image_dB );                                
colorbar; xlabel('Range(m)'); axis ij; ylabel('Doppler frequency (Hz)')                                           
title('Unfocused ISAR Image (no RA, no AF)'); axis xy; colormap('jet')
%% Range Alignment of Profiles 
% Range Align the HRR profiles using correlation method
[RA_HRR_profiles,shifts] = correlationRA(HRR_profiles);

% Plot range-aligned HRR Profiles
figure; imagesc(sb_HRR.G1.xaxis_downrange_m, 1:size(RA_HRR_profiles,1), ...
    20*log10(abs(RA_HRR_profiles))); 
xlabel('Range (m)');
ylabel('Profile Number');
title('Subset of HRR profiles (RA, no AF)');
colormap('jet');
colorbar;

% Plot Step Function
figure; plot(1:size(RA_HRR_profiles,1),shifts)
xlabel('Profile Number');
ylabel('Number of bin shifts')
title('Bin shifts per Range Profile');

%% Plot range-aligned ISAR image - no autofocus
% apply hamming window function
window = hamming(size(ProfilesToProcess,1));

ISAR_image = fftshift(fft(RA_HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension                    

ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 

% Plot RA ISAR image
figure; imagesc( Range_axis, DopplerAxis_Hz, ISAR_image_dB );                                
colorbar; xlabel('Range(m)'); axis ij; ylabel('Doppler frequency (Hz)')                                           
title('ISAR Image (RA, no AF)'); axis xy; colormap('jet')


