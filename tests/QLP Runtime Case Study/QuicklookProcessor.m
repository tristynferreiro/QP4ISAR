% ISAR quicklook processor
    % This was the initial version of the QLP used for testing the code
    % before porting it over to the APP.
%% 1 Get input parameters
CPI = input("What CPTWL would you like?\n");
overlap_percentage = input("What frame overlap % would you like?\n");

%% 2 Load CSIR dataset 
load('DAP_2010-10-14_09-43-33_010_zayaan_inbound_singlebin_P455_G1_sb_HRR.mat');
HRRProfiles_all = sb_HRR.G1.HRR_NoMC_calib.';
HRRProfiles_all = HRRProfiles_all(1:end,:); % adjust to data size wanted

profile_repetition_freq =  1/sb_HRR.G1.pattern_time; 
Range_axis = sb_HRR.G1.xaxis_downrange_m; 

[num_profiles,num_range_bins] = size(HRRProfiles_all);

% Check CPTWL is within range
while(CPI>num_profiles)
    CPI = input(sprintf("Error: Chosen CPTWL exceeds max value of %d. What CPTWL would you like?\n",num_range_bins));
end
%% 3 Slide window through all profiles (create frames)
ISAR_image=0;
% Step 1: Define parameters
overlap = overlap_percentage/100 * CPI;

% Step 2: Calculate the number of frames
num_frames = floor(1 + (num_profiles-CPI)/(CPI-overlap)); 

vidObj = VideoWriter ('ISARmovieZayaan1','MPEG-4');
vidObj.Quality = 100;
vidObj.FrameRate = 3;                   % Frames per second
open(vidObj); 

figure;

% Step 3: Get all unfocused frames
for frame = 1: num_frames
    startProfile = floor((frame-1) * (CPI-overlap) + 1 ); 
    stopProfile = floor(startProfile + CPI-1);
    ProfilesToProcess = stopProfile - startProfile;
    DopplerAxis_Hz = (-ProfilesToProcess/2:1:ProfilesToProcess/2-1)*profile_repetition_freq/ProfilesToProcess;

    %% 3.1 Unfocused HRRPs and ISAR image - no range alignment, no autofocus
    % Plot Unaligned HRR Profiles (single frame)
    HRR_profiles = circshift(HRRProfiles_all(startProfile:stopProfile, :), [0 50]);

    % apply hamming window function
    window = hamming(size(ProfilesToProcess,1));
    ISAR_image = fftshift(fft(HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension  
    % contrast value
    contrast_unfocused = imageContrast(ISAR_image);
    %fprintf("\nThe contrast IC value is %.4f",contrast_unfocused);
    %% 3.2 Range-aligned HRRPs and ISAR image - no autofocus
    % Range-align the HRR profiles
    ref_profile_number = 1; % arbitrarily selected
    RA_HRR_profiles = HaywoodRA(HRR_profiles,ref_profile_number);      
    %% 3.3 Autofocus Profiles
    % Autofocus the HRR profiles
    AF_RA_HRR_profiles = YuanAF(RA_HRR_profiles);

    % Plot focused ISAR image
    % apply hamming window function
    [num_profiles,num_range_bins] = size(RA_HRR_profiles);
    window = repmat(hamming(num_profiles),1,num_range_bins);

    ISAR_image = fftshift(fft(AF_RA_HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension                    

    ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 

    contrast_focused = roundn(imageContrast(ISAR_image),2);
    %fprintf("\nThe contrast IC value is %.4f",contrast_focused);

    %% Plot ISAR image
    cla; % clear the video axes
    imagesc( Range_axis, DopplerAxis_Hz, ISAR_image_dB );                                
    colorbar; xlabel('Range(m)'); axis ij; ylabel('Doppler frequency (Hz)')                                           
    title(['ISAR ' num2str(frame) ' of ' num2str(num_frames)  ' - C = ' num2str(contrast_focused)]); axis xy; colormap('jet');
    axis xy;
    pause(0.2);
    
    hold on;
    drawnow;
    ISARmov = getframe(gcf);
    writeVideo(vidObj, ISARmov); 
    hold off;
    
end
close(vidObj);
