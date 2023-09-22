% ISAR quicklook processor
%% 1 Get input parameters
CPTWL = input("What CPTWL would you like?\n");
overlap_percentage = input("What frame overlap % would you like?\n");
RA_selection = input("Range-alignment?\n 1 = Simple Correlation\n 2 = Haywood\n");
AF_selection = input("Autofocus?\n 1 = Yuan \n 2 = Haywood\n");

%% 2 Load CSIR dataset 
load('DAP_2010-10-14_09-43-33_010_zayaan_inbound_singlebin_P455_G1_sb_HRR.mat');
HRRProfiles_all = sb_HRR.G1.HRR_NoMC_calib.';
HRRProfiles_all = HRRProfiles_all(2400:end,:); % adjust to data size wanted

profile_repetition_freq =  1/sb_HRR.G1.pattern_time; 

[num_range_bins,num_profiles] = size(HRRProfiles_all);
Range_axis = sb_HRR.G1.xaxis_downrange_m; 

% Plot the HRR profiles (all)
% figure; imagesc(sb_HRR.G1.xaxis_downrange_m, 1:size(HRRProfilesAll,1), 20*log10(abs(HRRProfilesAll))); 
% xlabel('Range (m)');
% ylabel('Profile Number');
% title('All HRR profiles');
% colormap('jet');
% colorbar;

%% 3 Slide window through all profiles (create frames)
ISAR_image=0;
% Step 1: Define parameters
overlap = overlap_percentage/100 * CPTWL;

% Step 2: Calculate the number of frames
num_frames = floor(1 + (num_range_bins-CPTWL)/(CPTWL-overlap)); 

vidObj = VideoWriter ('ISARmovieRAHaywood.avi');
vidObj.FrameRate = 2;                   % Frames per second
open(vidObj); 

figure;

tic
% Step 3: Get all unfocused frames
for frame = 1: num_frames
    startProfile = floor((frame-1) * (CPTWL-overlap) + 1 ); 
    stopProfile = floor(startProfile + CPTWL-1);
    
    ProfilesToProcess = stopProfile - startProfile;
    DopplerAxis_Hz = (-ProfilesToProcess/2:1:ProfilesToProcess/2-1)*ProfileRepetitionFreq/ProfilesToProcess;

    %% 3.1 Unfocused HRRPs and ISAR image - no range alignment, no autofocus
    % Plot Unaligned HRR Profiles (single frame)
    HRR_profiles = circshift(HRRProfilesAll(startProfile:stopProfile, :), [0 50]);
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
    contrast_unfocused = roundn(imageContrast(ISAR_image),-2);
   
    %ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 
    
    %% 3.2 Range-aligned HRRPs and ISAR image - no autofocus
    % Range-align the HRR profiles using selected method
    ref_profile_number = 1;
    if(RA_selection == 1)
        [RA_HRR_profiles,shifts] = correlationRA(HRR_profiles,ref_profile_number);
        % % Plot Step Function
        % figure; plot(1:size(RA_HRR_profiles,1),shifts)
        % xlabel('Profile Number');
        % ylabel('Number of bin shifts')
        % title('Bin shifts per Range Profile');
    elseif (RA_selection == 2)
        [RA_HRR_profiles] = HaywoodRA(HRR_profiles,ref_profile_number);
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
    [num_profiles,num_range_bins] = size(HRR_profiles);
    window = repmat(hamming(num_profiles),1,num_range_bins);

    ISAR_image = fftshift(fft(RA_HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension                    

    % ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 

    % contrast = roundn(imageContrast(ISAR_image),2);
    % fprintf("\nThe contrast IC value is %.4f",contrast_beforeRA)

    %% 3.3 Autofocus Profiles
    % Range-align the HRR profiles using selected method
    if(AF_selection == 1)
        [AF_RA_HRR_profiles] = YuanAF(RA_HRR_profiles);
    elseif (AF_selection == 2)
        [AF_RA_HRR_profiles] = HaywoodAF(RA_HRR_profiles);

    end

    % Plot focused ISAR image
    % apply hamming window function
    window = repmat(hamming(num_profiles),1,num_range_bins);

    ISAR_image = fftshift(fft(AF_RA_HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension                    

    ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 

    contrast = roundn(imageContrast(ISAR_image),2);
    % fprintf("\nThe contrast IC value is %.4f",contrast_beforeRA)

    %% Plot ISAR image
    cla; % clear the video axes
    imagesc( Range_axis, DopplerAxis_Hz, ISAR_image_dB );                                
    colorbar; xlabel('Range(m)'); axis ij; ylabel('Doppler frequency (Hz)')                                           
    title(['ISAR ' num2str(frame) ' of ' num2str(num_frames)  ' - C = ' num2str(contrast)]); axis xy; colormap('jet')
    % pause(0.2);
    
    hold on;
    drawnow;
    ISARmov = getframe(gcf);
    writeVideo(vidObj, ISARmov); 
    hold off;
    
end
toc

close(vidObj);
