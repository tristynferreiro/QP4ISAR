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

% Plot the HRR profiles (all)
figure; imagesc(sb_HRR.G1.xaxis_downrange_m, 1:size(HRRProfilesAll,1), 20*log10(abs(HRRProfilesAll))); 
xlabel('Range (m)');
ylabel('Profile Number');
title('All HRR profiles');
colormap('jet');
colorbar;

%% Slide window through all profiles (create frames)
% Step 1: Define parameters
CPTWL = 400;
overlap = 0.8 * CPTWL;

% Step 2: Calculate the number of frames
num_frames = 1 + (NumRangeBins-CPTWL)/(CPTWL-overlap); 

% Step 3: Implement sliding window
for frame = 1: num_frames
    startProfile = (frame-1) * (CPTWL-overlap) + 1 ; 
    stopProfile = startProfile + CPTWL-1;
    
    ProfilesToProcess = stopProfile - startProfile;
    DopplerAxis_Hz = (-ProfilesToProcess/2:1:ProfilesToProcess/2-1)*ProfileRepetitionFreq/ProfilesToProcess;

    % Plot Unaligned HRR Profiles (single frame)
    HRR_profiles = circshift(HRRProfilesAll(startProfile:stopProfile, :), [0 50]);
    figure; imagesc(sb_HRR.G1.xaxis_downrange_m, 1:size(HRR_profiles,1), 20*log10(abs(HRR_profiles))); 
    xlabel('Range (m)');
    ylabel('Profile Number');
    title('Subset of HRR profiles');
    colormap('jet');
    colorbar;


    % %% Plot unfocused ISAR image - no range alignment, no autofocus
    % % apply hamming window function
    % window = hamming(size(ProfilesToProcess,1));
    % 
    % ISAR_image = fftshift(fft(HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension                    
    % 
    % ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 
    % 
    % % Plot unfocused ISAR image
    % figure; imagesc( Range_axis, DopplerAxis_Hz, ISAR_image_dB );                                
    % colorbar; xlabel('Range(m)'); axis ij; ylabel('Doppler frequency (Hz)')                                           
    % title('Unfocused ISAR Image (no RA, no AF)'); axis xy; colormap('jet')
   
end





