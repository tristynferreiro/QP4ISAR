function [contrast_unfocused,contrast_focused] = QLPDesignImageContrastTesting(CPTWL,frame_center,RA_selection,AF_selection)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% 2 Load CSIR dataset 
load('DAP_2010-10-14_09-43-33_010_zayaan_inbound_singlebin_P455_G1_sb_HRR.mat');

HRR_profiles_all = sb_HRR.G1.HRR_NoMC_calib.';       % All HRR profiles
Profile_Repetition_Freq =  1/sb_HRR.G1.pattern_time; % PRF
Range_axis = sb_HRR.G1.xaxis_downrange_m; 

%% 3 Choose middle profile
StartProfile = frame_center - CPTWL/2;                           
StopProfile = frame_center+ CPTWL/2 - 1;
ProfilesToProcess = StopProfile - StartProfile;
DopplerAxis_Hz = (-ProfilesToProcess/2:1:ProfilesToProcess/2-1)*Profile_Repetition_Freq/ProfilesToProcess;
%% 4 Unfocused Image
% Unaligned HRR Profiles
HRR_profiles = circshift(HRR_profiles_all(StartProfile:StopProfile, :), [0 50 ]);

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
ref_profile_number =1; % arbitrarily selected
if(RA_selection == 1)
    RA_HRR_profiles = correlationRA(HRR_profiles,ref_profile_number);
elseif (RA_selection == 2)
    RA_HRR_profiles = HaywoodRA(HRR_profiles,ref_profile_number);
end

%% Autofocus of Profiles
if(AF_selection ~= 0)
    % Autofocus the HRR profiles using selected method
    if(AF_selection == 1)
        AF_RA_HRR_profiles = YuanAF(RA_HRR_profiles);
    elseif (AF_selection == 2)
        AF_RA_HRR_profiles = HaywoodAF(RA_HRR_profiles);
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
    contrast_focused = imageContrast(ISAR_image);
    fprintf("\nThe contrast IC value after AF is %.4f",contrast_focused);
end

end