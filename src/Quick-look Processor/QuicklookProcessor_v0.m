% ISAR quicklook processor
    % This was the initial version of the QLP used for during development 
    % and still contains many intemediary steps and plots used.
clear all;
%% 1 Get input parameters
%CPTWL = input("What CPTWL would you like?\n");
CPI = 128;
overlap_percentage = 20;
%overlap_percentage = input("What frame overlap % would you like?\n");
RA_selection = 2;
AF_selection = 1;
% RA_selection = input("Range-alignment?\n 0= none\n 1 = Simple Correlation\n 2 = Haywood\n");
% AF_selection = input("Autofocus?\n 0=none\n 1 = Yuan \n 2 = Haywood\n");

tBegin = tic;
%% 2 Load CSIR dataset 
load('DAP_2010-10-14_09-43-33_010_zayaan_inbound_singlebin_P455_G1_sb_HRR.mat');
% load('DAP_2010-10-14_18-59-55_014_Zayaan_P870_G1_sb_HRR.mat'); % Zayaan 2
%load('DAP_2010-10-12_11-09-28_004_Tugboat_P874_G1_sb_HRR.mat'); % Tugboat
HRRProfiles_all = sb_HRR.G1.HRR_NoMC_calib.';
HRRProfiles_all = HRRProfiles_all(1:end,:); %2400 % adjust to data size wanted

profile_repetition_freq =  1/sb_HRR.G1.pattern_time; 
profile_repetition_time = sb_HRR.G1.pattern_time;
Range_axis = sb_HRR.G1.xaxis_downrange_m; 

[num_profiles,num_range_bins] = size(HRRProfiles_all);

% Plot the HRR profiles (all)
% figure; imagesc(num_range_bins, 1:size(HRRProfilesAll,1), 20*log10(abs(HRRProfilesAll))); 
% xlabel('Range (m)');
% ylabel('Profile Number');
% title('All HRR profiles');
% colormap('jet');
% colorbar;

% Calculate CPTWL
acquisition_time = profile_repetition_time * num_profiles;
fprintf("The data measurement window is %.4f (s).\n",acquisition_time)

%% 3 Slide window through all profiles (create frames)
ISAR_image=0;

% Step 1: Define parameters
overlap = overlap_percentage/100 * CPI;

% Step 2: Calculate the number of frames
num_frames = floor(1 + (num_profiles-CPI)/(CPI-overlap)); 

vidObj = VideoWriter ('ISARmovie.avi');
vidObj.FrameRate = 2;                   % Frames per second
open(vidObj); 

figure;

% Step 3: Get all unfocused frames
for frame = 1: num_frames
    %fprintf("Frame is %d\n",frame);

    startProfile = floor((frame-1) * (CPI-overlap) + 1 ); 
    stopProfile = floor(startProfile + CPI-1);
    
    ProfilesToProcess = (stopProfile - startProfile)+1;
    DopplerAxis_Hz = (-ProfilesToProcess/2:1:ProfilesToProcess/2-1)*profile_repetition_freq/ProfilesToProcess;

    %% 3.1 Unfocused HRRPs and ISAR image - no range alignment, no autofocus
    % Plot Unaligned HRR Profiles (single frame)
    HRR_profiles = circshift(HRRProfiles_all(startProfile:stopProfile, :), [0 50]);
    % figure; imagesc(sb_HRR.G1.xaxis_downrange_m, 1:size(HRR_profiles,1), 20*log10(abs(HRR_profiles))); 
    % xlabel('Range (m)');
    % ylabel('Profile Number');
    % title('Subset of HRR profiles');
    % colormap('jet');
    % colorbar;

    % apply hamming window function
    window = hamming(size(ProfilesToProcess,1));
    ISAR_image = fftshift(fft(HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension  
    % contrast value
    % contrast_unfocused = imageContrast(ISAR_image_dB);
    % fprintf("\nThe contrast IC value of unfocused image is %.4f",contrast_unfocused)
   
    %ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 
    
    %% 3.2 Range-aligned HRRPs and ISAR image - no autofocus
    if(RA_selection~=0)
        % Range-align the HRR profiles using selected algorithm
        ref_profile_number = 1;
        if(RA_selection == 1)
            RA_HRR_profiles = correlationRA(HRR_profiles,ref_profile_number);
        elseif (RA_selection == 2)
            RA_HRR_profiles = HaywoodRA(HRR_profiles,ref_profile_number);
        end
    
        % % Plot range-aligned HRR Profiles
        % figure; imagesc(sb_HRR.G1.xaxis_downrange_m, 1:size(RA_HRR_profiles,1), ...
        %     20*log10(abs(RA_HRR_profiles))); 
        % xlabel('Range (m)');
        % ylabel('Profile Number');
        % title('Subset of HRR profiles (RA, no AF)');
        % colormap('jet');
        % colorbar;
    
        % Plot range-aligned ISAR image - no autofocus
        % apply hamming window function
        [num_range_bins,num_profiles] = size(RA_HRR_profiles);
        window = repmat(hamming(num_range_bins),1,num_profiles);
    
        ISAR_image = fftshift(fft(RA_HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension                    
    
        ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 
    
        % contrast_afterRA = imageContrast(ISAR_image);
        % fprintf("\nThe contrast IC value after RA is %.4f",contrast_afterRA)
    end
    %% 3.3 Autofocus Profiles
    if(AF_selection~=0)
        % Autofocus the HRR profiles using selected method
        if(AF_selection == 1)
            AF_RA_HRR_profiles = YuanAF_v2(RA_HRR_profiles);
        elseif (AF_selection == 2)
            AF_RA_HRR_profiles = HaywoodAF(RA_HRR_profiles);
    
        end
    
        % Plot focused ISAR image
        % apply hamming window function
        window = repmat(hamming(num_range_bins),1,num_profiles);
    
        ISAR_image = fftshift(fft(AF_RA_HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension                    
    
        ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 
    
        contrast_afterAF = imageContrast(ISAR_image);
        fprintf("\nThe contrast IC value after AF is %.4f",contrast_afterAF)
    end
    %% Plot ISAR image
    cla; % clear the video axes
    imagesc( Range_axis, DopplerAxis_Hz, ISAR_image_dB );                                
    colorbar; xlabel('Range(m)'); axis ij; ylabel('Doppler frequency (Hz)')                                           
    title(['ISAR ' num2str(frame) ' of ' num2str(num_frames)  ' - C = ' num2str(contrast_afterAF)]); axis xy; colormap('jet')
    axis xy;
    % pause(0.2);

    hold on;
    drawnow;
    ISARmov = getframe(gcf);
    writeVideo(vidObj, ISARmov); 
    hold off;
    
end
tEnd = toc;

close(vidObj);
