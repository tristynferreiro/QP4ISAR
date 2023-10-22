function [] = QLPDesignTimeTesting(CPTWL,overlap_percentage,RA_selection,AF_selection,video_name)
% Setup used for time testing for QLP design considerations.
%   This function was used to run time testing of all RA and AF
%   combinations considered for use in the final QLP design. 
%
%   Note: only the core functionality is included in this testing function,
%   i.e. no intemediary plots or image contrast calculations are included.

%% 1 Load CSIR dataset 
load('DAP_2010-10-14_09-43-33_010_zayaan_inbound_singlebin_P455_G1_sb_HRR.mat');
HRRProfiles_all = sb_HRR.G1.HRR_NoMC_calib.';
HRRProfiles_all = HRRProfiles_all(1:end,:); %2400 % adjust to data size wanted

profile_repetition_freq =  1/sb_HRR.G1.pattern_time; 
profile_repetition_time = sb_HRR.G1.pattern_time;
Range_axis = sb_HRR.G1.xaxis_downrange_m; 

[num_profiles,num_range_bins] = size(HRRProfiles_all);

%% 3 Slide window through all profiles (create frames)
ISAR_image = 0;

vidObj = VideoWriter (video_name);
vidObj.FrameRate = 2;                   % Frames per second
open(vidObj); 

figure;

% Step 1: Define parameters
overlap = overlap_percentage/100 * CPTWL;

% Step 2: Calculate the number of frames
num_frames = floor(1 + (num_profiles-CPTWL)/(CPTWL-overlap)); 

% Step 3: Get all unfocused frames
for frame = 1: num_frames
    %fprintf("Frame is %d\n",frame);

    startProfile = floor((frame-1) * (CPTWL-overlap) + 1 ); 
    stopProfile = floor(startProfile + CPTWL-1);
    
    ProfilesToProcess = (stopProfile - startProfile)+1;
    DopplerAxis_Hz = (-ProfilesToProcess/2:1:ProfilesToProcess/2-1)*profile_repetition_freq/ProfilesToProcess;

    %% 3.1 Unfocused HRRPs and ISAR image - no range alignment, no autofocus
    % Plot Unaligned HRR Profiles (single frame)
    HRR_profiles = circshift(HRRProfiles_all(startProfile:stopProfile, :), [0 50]);

    % apply hamming window function
    window = hamming(size(ProfilesToProcess,1));
    ISAR_image = fftshift(fft(HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension  
    
    %% 3.2 Range-aligned HRRPs and ISAR image - no autofocus
    if(RA_selection~=0)
        % Range-align the HRR profiles using selected algorithm
        ref_profile_number = 1;
        if(RA_selection == 1)
            RA_HRR_profiles = correlationRA(HRR_profiles,ref_profile_number);
        elseif (RA_selection == 2)
            RA_HRR_profiles = HaywoodRA(HRR_profiles,ref_profile_number);
        end
    
        % Plot range-aligned ISAR image - no autofocus
        % apply hamming window function
        [num_range_bins,num_profiles] = size(RA_HRR_profiles);
        window = repmat(hamming(num_range_bins),1,num_profiles);
    
        ISAR_image = fftshift(fft(RA_HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension                    
    
        ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 
    end
    %% 3.3 Autofocus Profiles
    if(AF_selection~=0)
        % Autofocus the HRR profiles using selected method
        if(AF_selection == 1)
            AF_RA_HRR_profiles = YuanAF(RA_HRR_profiles);
        elseif (AF_selection == 2)
            AF_RA_HRR_profiles = HaywoodAF(RA_HRR_profiles);
    
        end
    
        % Plot focused ISAR image
        % apply hamming window function
        window = repmat(hamming(num_range_bins),1,num_profiles);
    
        ISAR_image = fftshift(fft(AF_RA_HRR_profiles.*window,[],1),1); % Apply FFT in the slow-time dimension                    
    
        ISAR_image_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_image, 35);    % Limit dynamic range of ISAR image 
        
    end
    %% Plot Focused ISAR image
    contrast_after = imageContrast(ISAR_image); % calculate IC
    cla; % clear the video axes
    imagesc( Range_axis, DopplerAxis_Hz, ISAR_image_dB );                                
    colorbar; xlabel('Range(m)'); axis ij; ylabel('Doppler frequency (Hz)')                                           
    title(['ISAR ' num2str(frame) ' of ' num2str(num_frames)  ' - C = ' num2str(contrast_after)]); 
    axis xy; colormap('jet');
    pause(0.2);

    hold on;
    drawnow;
    ISARmov = getframe(gcf);
    writeVideo(vidObj, ISARmov); 
    hold off;
    
en
close(vidObj);
end