%close all; %clear all; clc; 

%% 1 Get input parameters
middle_profile = input("Middle Profile?\n 1 = 2464\n 2 = 4189\n 3 = 2970\n 4 = 3827\n");
RA_selection = input("Range-alignment?\n 0 = none\n 1 = Simple Correlation\n 2 = Haywood\n");
AF_selection = input("Autofocus?\n 0 = none\n 1 = Yuan \n 2 = Haywood\n");
%% 2 Load CSIR dataset 

load('DAP_2010-10-14_09-43-33_010_zayaan_inbound_singlebin_P455_G1_sb_HRR.mat');
HRRProfilesAll = sb_HRR.G1.HRR_NoMC_calib.';
ProfileRepetitionFreq =  1/sb_HRR.G1.pattern_time; 
NumRangeBins = size(HRRProfilesAll,1);
NumOfProfiles = size(HRRProfilesAll,2);
Range_axis = sb_HRR.G1.xaxis_downrange_m; 

% Plot the HRR profiles
% figure; imagesc(sb_HRR.G1.xaxis_downrange_m, 1:size(HRRProfilesAll,1), 20*log10(abs(HRRProfilesAll))); 
% xlabel('Range (m)');
% ylabel('Profile Number');
% title('All HRR profiles');
% colormap('jet');
% colorbar;
%% 3 Choose middle profile
% Select a subset of profiles
if middle_profile==1
    MiddleProfile = 2464;                                              % 2464, 4189, 2970, 3827
elseif middle_profile==2
    MiddleProfile = 4189;
elseif middle_profile==3
    MiddleProfile = 2970;
elseif middle_profile==4
    MiddleProfile = 3827;
end

CPTWL = 128;   
StartProfile = MiddleProfile - CPTWL/2;                           
StopProfile = MiddleProfile + CPTWL/2 - 1;
ProfilesToProcess = StopProfile - StartProfile;
DopplerAxis_Hz = (-ProfilesToProcess/2:1:ProfilesToProcess/2-1)*ProfileRepetitionFreq/ProfilesToProcess;
%% 4 Unfocused Image

% Plot Unaligned HRR Profiles
HRR_profiles = circshift(HRRProfilesAll(StartProfile:StopProfile, :), [0 50 ]);
figure; imagesc(sb_HRR.G1.xaxis_downrange_m, 1:size(HRR_profiles,1), 20*log10(abs(HRR_profiles))); 
xlabel('Range (m)');
ylabel('Profile Number');
title('Subset of HRR profiles');
colormap('jet');
colorbar;
% matlab2tikz() % Save the figure as LaTeX compatible plot

% Unfocused ISAR image - no range alignment, no autofocus
% apply hamming window function
WindowMatrix = repmat(size(ProfilesToProcess,1),1, size(ProfilesToProcess,2));
ISAR_image = fftshift(fft(HRR_profiles.*WindowMatrix,[],1),1); % Apply FFT in the slow-time dimension                    
ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 

% Plot unfocused ISAR image
figure; imagesc( Range_axis, DopplerAxis_Hz, ISAR_image_dB );                                
colorbar; xlabel('Range(m)'); axis ij; ylabel('Doppler frequency (Hz)')                                           
title('Unfocused ISAR Image (no RA, no AF)'); axis xy; colormap('jet')

% matlab2tikz() % Save the figure as LaTeX compatible plot
%% 5 Range Alignment of Profiles 
if RA_selection ~=0
    % Range Align the HRR profiles using user selected algorithm
    ref_profile_number =1;
    if(RA_selection == 1)
        [RA_HRR_profiles,shifts] = correlationRA(HRR_profiles,ref_profile_number);
        % Plot Step Function
        % figure; plot(1:size(RA_HRR_profiles,1),shifts)
        % xlabel('Profile Number');
        % ylabel('Number of bin shifts')
        % title('Bin shifts per Range Profile');
        % matlab2tikz() % Save the figure as LaTeX compatible plot
    elseif (RA_selection == 2)
        [RA_HRR_profiles] = HaywoodRA(HRR_profiles,ref_profile_number);
    end

    
    % Plot range-aligned HRR Profiles
    figure; imagesc(sb_HRR.G1.xaxis_downrange_m, 1:size(RA_HRR_profiles,1), ...
        20*log10(abs(RA_HRR_profiles))); 
    xlabel('Range (m)');
    ylabel('Profile Number');
    title('Range-aligned HRR profiles');
    colorbar;
    colormap('jet');
    % matlab2tikz() % Save the figure as LaTeX compatible plot
    
    % Plot range-aligned ISAR image - no autofocus
    WindowMatrix = repmat(size(ProfilesToProcess,1),1, size(ProfilesToProcess,2)); % apply hamming window function
    ISAR_image = fftshift(fft(RA_HRR_profiles.*WindowMatrix,[],1),1); % Apply FFT in the slow-time dimension                    
    ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 
    
    % Plot RA ISAR image
    figure; imagesc(Range_axis, DopplerAxis_Hz, ISAR_image_dB);                                
    xlabel('Range(m)'); 
    ylabel('Doppler frequency (Hz)')                                           
    title('Range-aligned ISAR image');
    axis xy; 
    colorbar; 
    colormap('jet')
    % matlab2tikz() % Save the figure as LaTeX compatible plot
    
    % Calculate contrast value
    % contrast_beforeRA = imageContrast(ISAR_image);
    % fprintf("\nThe contrast IC value is %.4f",contrast_beforeRA)
end

%% Autofocus of Profiles
if(AF_selection ~= 0)
    % Autofocus the HRR profiles using selected method
    if(AF_selection == 1)
        [AF_RA_HRR_profiles] = YuanAF(RA_HRR_profiles);
    elseif (AF_selection == 2)
        [AF_RA_HRR_profiles] = HaywoodAF(RA_HRR_profiles);
    end
    
    % Plot ISAR image
    WindowMatrix = repmat(size(ProfilesToProcess,1),1, size(ProfilesToProcess,2));
    ISAR_linear = fftshift(fft(AF_RA_HRR_profiles.*WindowMatrix, [], 1),1);
    ISAR_linear_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_linear,35);
    
    figure;
    imagesc(RangeAxis, DopplerAxis_Hz, ISAR_linear_dB);
    xlabel('Range (m)');
    ylabel('Doppler frquency (Hz)');
    title('RA and AF Focused ISAR image');
    colorbar;
    colormap('jet');
    axis xy;
    % matlab2tikz() % Save the figure as LaTeX compatible plot
end
