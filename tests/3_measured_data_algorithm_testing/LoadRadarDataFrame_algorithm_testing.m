% This is a testing script that was used during validation of the RA and AF
    % algorithms. 
    % 
    % It loads in a MATLAB structure containg all HRR profiles. Based on
    % user selection, it greates a single data frame (group of HRR profiles
    % and range bins). Again, based on user selection, RA and AF algorithms
    % are performed on the data frame. 
    % 
    % Inputs:
    %   - MATLAB data structure containg HRR profiles
    %   - User data frame selection: which frame of data to use
    %   - User RA selection: Which RA algorithm to use
    %   - User AF selection: Which AF algorithm to use
    %
    % Outputs:
    %   - HRR profiles of the selected frame.
    %   - ISAR images before and after each algorithm is applied.
    
close all; clear all; clc; 

%% 1 Get input parameters
frame_center = input("Frame Center?\n 1 = 2464\n 2 = 4189\n 3 = 2970\n 4 = 3827\n");
RA_selection = input("Range-alignment?\n 0 = none\n 1 = Simple Correlation\n 2 = Haywood\n");
AF_selection = input("Autofocus?\n 0 = none\n 1 = Yuan \n 2 = Haywood\n");
%% 2 Load CSIR dataset 

load('DAP_2010-10-14_09-43-33_010_zayaan_inbound_singlebin_P455_G1_sb_HRR.mat');
HRR_profiles_all = sb_HRR.G1.HRR_NoMC_calib.';       % All HRR profiles
Profile_Repetition_Freq =  1/sb_HRR.G1.pattern_time; % PRF
Num_range_bins = size(HRR_profiles_all,1);
Num_profiles = size(HRR_profiles_all,2);
Range_axis = sb_HRR.G1.xaxis_downrange_m; 

% Plot the HRR profiles
% figure; imagesc(sb_HRR.G1.xaxis_downrange_m, 1:size(HRRProfilesAll,1), 20*log10(abs(HRRProfilesAll))); 
% xlabel('Range (m)'); ylabel('Profile Number');
% title('All HRR profiles');
% colormap('jet'); colorbar;

%% 3 Choose middle profile
% Select a subset of profiles
if frame_center==1
    Middle_of_frame = 2464;                                              % 2464, 4189, 2970, 3827
elseif frame_center==2
    Middle_of_frame = 4189;
elseif frame_center==3
    Middle_of_frame = 2970;
elseif frame_center==4
    Middle_of_frame = 3827;
end

CPI = 128;   
StartProfile = Middle_of_frame - CPI/2;                           
StopProfile = Middle_of_frame + CPI/2 - 1;
ProfilesToProcess = StopProfile - StartProfile;
DopplerAxis_Hz = (-ProfilesToProcess/2:1:ProfilesToProcess/2-1)*Profile_Repetition_Freq/ProfilesToProcess;
%% 4 Unfocused Image

% Plot Unaligned HRR Profiles
HRR_profiles = circshift(HRR_profiles_all(StartProfile:StopProfile, :), [0 50 ]);
% figure; imagesc(1:size(HRR_profiles,1), 1:size(HRR_profiles,2), 20*log10(abs(HRR_profiles))); 
% xlabel('Range Bin Number'); ylabel('Profile Number');
% title('Subset of HRR profiles');
% colormap('jet'); %colorbar;
% matlab2tikz() % Save the figure as LaTeX compatible plot

% Unfocused ISAR image - no range alignment, no autofocus
% apply hamming window function
window_matrix = repmat(size(ProfilesToProcess,1),1, size(ProfilesToProcess,2));
ISAR_image = fftshift(fft(HRR_profiles.*window_matrix,[],1),1); % Apply FFT in the slow-time dimension                    
ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 

% Plot unfocused ISAR image
figure; imagesc( Range_axis, DopplerAxis_Hz, ISAR_image_dB );                                
xlabel('Range(m)'); ylabel('Doppler frequency (Hz)')                                           
title('Unfocused ISAR Image (no RA, no AF)'); 
axis xy; colormap('jet'); %colorbar;
% matlab2tikz() % Save the figure as LaTeX compatible plot

% Calculate contrast value
contrast_unfocused = imageContrast(ISAR_image);
fprintf("\nThe contrast IC value of the unfocused image is %.4f",contrast_unfocused);
%% 5 Range Alignment of Profiles 
if RA_selection ~=0
    % Range Align the HRR profiles using user selected algorithm
    ref_profile_number =1;
    if(RA_selection == 1)
        RA_HRR_profiles = correlationRA(HRR_profiles,ref_profile_number);
    elseif (RA_selection == 2)
        RA_HRR_profiles = HaywoodRA(HRR_profiles,ref_profile_number);
    end
    
    % Plot range-aligned HRR Profiles
    % figure; imagesc(1:size(HRR_profiles,1), 1:size(HRR_profiles,2), ...
    %     20*log10(abs(RA_HRR_profiles))); 
    % xlabel('Range Bin Number'); ylabel('Profile Number');
    % title('Range-aligned HRR profiles');
    % colorbar; colormap('jet');
    % matlab2tikz() % Save the figure as LaTeX compatible plot
    
    % Plot range-aligned ISAR image - no autofocus
    window_matrix = repmat(size(ProfilesToProcess,1),1, size(ProfilesToProcess,2)); % apply hamming window function
    ISAR_image = fftshift(fft(RA_HRR_profiles.*window_matrix,[],1),1); % Apply FFT in the slow-time dimension                    
    ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 
    
    % Plot RA ISAR image
    figure; imagesc(Range_axis, DopplerAxis_Hz, ISAR_image_dB);                                
    xlabel('Range(m)'); ylabel('Doppler frequency (Hz)')                                           
    title('Range-aligned ISAR image');
    axis xy; colorbar; colormap('jet')
    % matlab2tikz() % Save the figure as LaTeX compatible plot
    
    % Calculate contrast value
    contrast_afterRA = imageContrast(ISAR_image);
    fprintf("\nThe contrast IC value after RA is %.4f",contrast_afterRA);
end

%% Autofocus of Profiles
if(AF_selection ~= 0)
    % Autofocus the HRR profiles using selected method
    autofocus = "";
    if(AF_selection == 1)
        AF_RA_HRR_profiles = YuanAF(RA_HRR_profiles);
        autofocus = "Yuan";
    elseif (AF_selection == 2)
        AF_RA_HRR_profiles = HaywoodAF(RA_HRR_profiles);
        autofocus = "Haywood";
    end
    
    % Plot ISAR image
    window_matrix = repmat(size(ProfilesToProcess,1),1, size(ProfilesToProcess,2));
    ISAR_image = fftshift(fft(AF_RA_HRR_profiles.*window_matrix, [], 1),1);
    ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image,35);
    
    figure; imagesc(Range_axis, DopplerAxis_Hz, ISAR_image_dB);
    xlabel('Range (m)'); ylabel('Doppler frquency (Hz)');
    title('RA and AF Focused ISAR image');
    colorbar; colormap('jet'); axis xy;
    % matlab2tikz() % Save the figure as LaTeX compatible plot

    % Calculate contrast value
    contrast_afterAF = imageContrast(ISAR_image);
    fprintf("\nThe contrast IC value after AF is %.4f",contrast_afterAF);
end

