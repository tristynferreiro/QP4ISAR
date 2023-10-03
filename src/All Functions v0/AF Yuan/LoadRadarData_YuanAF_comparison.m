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

% % Plot the HRR profiles
% figure; imagesc(sb_HRR.G1.xaxis_downrange_m, 1:size(HRRProfilesAll,1), 20*log10(abs(HRRProfilesAll))); 
% xlabel('Range (m)');
% ylabel('Profile Number');
% title('All HRR profiles');
% colormap('jet');
% colorbar;

% Select a subset of profiles
MiddleProfile = 2464;                                              % 2464, 4189, 2970, 3827
CPTWL = 128;   
StartProfile = MiddleProfile - CPTWL/2;                           
StopProfile = MiddleProfile + CPTWL/2 - 1;
ProfilesToProcess = StopProfile - StartProfile;
DopplerAxis_Hz = (-ProfilesToProcess/2:1:ProfilesToProcess/2-1)*ProfileRepetitionFreq/ProfilesToProcess;

% % Plot Unaligned HRR Profiles
HRR_profiles = circshift(HRRProfilesAll(StartProfile:StopProfile, :), [0 50 ]);
% figure; imagesc(sb_HRR.G1.xaxis_downrange_m, 1:size(HRR_profiles,1), 20*log10(abs(HRR_profiles))); 
% xlabel('Range (m)');
% ylabel('Profile Number');
% title('Subset of HRR profiles');
% colormap('jet');
% colorbar;

%% Plot unfocused ISAR image - no range alignment, no autofocus
% apply hamming window function
window = hamming(size(ProfilesToProcess,1));

ISAR_image = fftshift(fft(HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension                    

ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 

% Plot unfocused ISAR image
figure; imagesc( Range_axis, DopplerAxis_Hz, ISAR_image_dB );                                
colorbar; xlabel('Range(m)'); axis ij; ylabel('Doppler frequency (Hz)')                                           
title('Unfocused ISAR Image (no RA, no AF)'); axis xy; colormap('jet')
%% Simple correlation Range Alignment of Profiles 
% Range Align the HRR profiles using correlation method
ref_profile_number = 1;
[corrRA_HRR_profiles,shifts] = correlationRA(HRR_profiles,ref_profile_number);

% Plot simple correlation range-aligned ISAR image - no autofocus
% apply hamming window function
[numProfiles,numRangeBins] = size(corrRA_HRR_profiles);
window = repmat(hamming(numProfiles),1,numRangeBins);

ISAR_image_corrRA = fftshift(fft(corrRA_HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension                    

ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image_corrRA, 35);    % Limit dynamic range of ISAR image 

% Plot RA ISAR image
figure; imagesc( Range_axis, DopplerAxis_Hz, ISAR_image_dB );                                
colorbar; xlabel('Range(m)'); axis ij; ylabel('Doppler frequency (Hz)')                                           
title('ISAR Image (corr RA, no AF)'); axis xy; colormap('jet')
%% Contrast comparison
contrast_beforeRA = imageContrast(ISAR_image);
contrast_afterRA = imageContrast(ISAR_image_corrRA);
% Calculate focus improvement factor
contrast_improvement = contrast_afterRA/contrast_beforeRA;
fprintf("\nThe contrast improvement factor after simple corr RA is %.4f",contrast_improvement)
%% Haywood Range Alignment of Profiles 
% Range Align the HRR profiles using Haywood method
ref_profile_number = 1;
HaywoodRA_HRR_profiles = HaywoodRA(HRR_profiles, ref_profile_number);

% Plot Haywood range-aligned ISAR image - no autofocus
% apply hamming window function
window = repmat(hamming(numProfiles),1,numRangeBins);

ISAR_image_HaywoodRA = fftshift(fft(HaywoodRA_HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension                    

ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image_HaywoodRA, 35); % Limit dynamic range of ISAR image 

% Plot RA ISAR image
figure; imagesc( Range_axis, DopplerAxis_Hz, ISAR_image_dB );                                
colorbar; xlabel('Range(m)'); axis ij; ylabel('Doppler frequency (Hz)')                                           
title('ISAR Image (Haywood RA, no AF)'); axis xy; colormap('jet')
%% Contrast comparison
contrast_beforeRA = imageContrast(ISAR_image);
contrast_afterRA = imageContrast(ISAR_image_HaywoodRA);
% Calculate focus improvement factor
contrast_improvement = contrast_afterRA/contrast_beforeRA;
fprintf("\nThe contrast improvement factor after Haywood RA is %.4f",contrast_improvement)
%% Autofocus of Profiles
% Apply Yuan autofocus to the RA HRR profiles using with noise
noise_AF_corrRA_HRR_profiles = YuanAF_noise(corrRA_HRR_profiles);
noise_AF_HaywoodRA_HRR_profiles = YuanAF_noise(HaywoodRA_HRR_profiles);

% Apply Yuan autofocus to the RA HRR profiles using
AF_corrRA_HRR_profiles = YuanAF(corrRA_HRR_profiles);
AF_HaywoodRA_HRR_profiles = YuanAF(HaywoodRA_HRR_profiles);

%% Plot focused ISAR image - with noise
% apply hamming window function
window = repmat(hamming(numProfiles),1,numRangeBins);

noise_focused_corrRA_ISAR_image = fftshift(fft(noise_AF_corrRA_HRR_profiles .*window,[],1),1); % Apply FFT in the slow-time dimension                    

focused_ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(noise_focused_corrRA_ISAR_image, 35); % Limit dynamic range of ISAR image 

% Plot AF (corr RA) ISAR image
figure; imagesc( Range_axis, DopplerAxis_Hz, focused_ISAR_image_dB );                                
colorbar; xlabel('Range(m)'); axis ij; ylabel('Doppler frequency (Hz)')                                           
title('ISAR Image (corr RA, noise AF)'); axis xy; colormap('jet')

% apply hamming window function
window = repmat(hamming(numProfiles),1,numRangeBins);

noise_focused_HaywoodRA_ISAR_image = fftshift(fft(noise_AF_HaywoodRA_HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension                    

focused_ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(noise_focused_HaywoodRA_ISAR_image, 35); % Limit dynamic range of ISAR image 

% Plot AF (Haywood RA) ISAR image
figure; imagesc( Range_axis, DopplerAxis_Hz, focused_ISAR_image_dB );                                
colorbar; xlabel('Range(m)'); axis ij; ylabel('Doppler frequency (Hz)')                                           
title('ISAR Image (Haywood RA, noise AF)'); axis xy; colormap('jet')
% apply hamming window function
window = repmat(hamming(numProfiles),1,numRangeBins);

focused_corrRA_ISAR_image = fftshift(fft(AF_corrRA_HRR_profiles .*window,[],1),1); % Apply FFT in the slow-time dimension                    

focused_ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(focused_corrRA_ISAR_image, 35); % Limit dynamic range of ISAR image 

% Plot AF (corr RA) ISAR image
figure; imagesc( Range_axis, DopplerAxis_Hz, focused_ISAR_image_dB );                                
colorbar; xlabel('Range(m)'); axis ij; ylabel('Doppler frequency (Hz)')                                           
title('ISAR Image (corr RA, AF)'); axis xy; colormap('jet')

% apply hamming window function
window = repmat(hamming(numProfiles),1,numRangeBins);

focused_HaywoodRA_ISAR_image = fftshift(fft(AF_HaywoodRA_HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension                    

focused_ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(focused_HaywoodRA_ISAR_image, 35); % Limit dynamic range of ISAR image 

% Plot AF (Haywood RA) ISAR image
figure; imagesc( Range_axis, DopplerAxis_Hz, focused_ISAR_image_dB );                                
colorbar; xlabel('Range(m)'); axis ij; ylabel('Doppler frequency (Hz)')                                           
title('ISAR Image (Haywood RA, AF)'); axis xy; colormap('jet')
%% Contrast comparison - with noise
contrast_before = imageContrast(ISAR_image);
contrast_afterAF = imageContrast(noise_focused_corrRA_ISAR_image);
% Calculate focus improvement factor
contrast_improvement = contrast_afterAF/contrast_before;
fprintf("\nThe contrast improvement factor after noise Yuan AF (corr RA) is %.4f",contrast_improvement)

contrast_before = imageContrast(ISAR_image);
contrast_afterAF = imageContrast(noise_focused_HaywoodRA_ISAR_image);
% Calculate focus improvement factor
contrast_improvement = contrast_afterAF/contrast_before;
fprintf("\nThe contrast improvement factor after noise Yuan AF (Haywood RA) is %.4f",contrast_improvement)

% Without noise
contrast_before = imageContrast(ISAR_image);
contrast_afterAF = imageContrast(focused_corrRA_ISAR_image);
% Calculate focus improvement factor
contrast_improvement = contrast_afterAF/contrast_before;
fprintf("\nThe contrast improvement factor after Yuan AF (corr RA) is %.4f",contrast_improvement)

contrast_before = imageContrast(ISAR_image);
contrast_afterAF = imageContrast(focused_HaywoodRA_ISAR_image);
% Calculate focus improvement factor
contrast_improvement = contrast_afterAF/contrast_before;
fprintf("\nThe contrast improvement factor after Yuan AF (Haywood RA) is %.4f",contrast_improvement)

