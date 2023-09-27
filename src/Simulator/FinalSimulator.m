
% Version 3 - include target translation and rotational motion 
% Revision: allows for adjustment of the scatterer voltages

clear all; close all; clc;
%% 1 User Input Parameters
target_velocity = input("What is the target's translational velocity (m/s)?\n");
rotation_angle = input("What is the rotation angle (degrees)?\n");
RA_selection = input("Range-alignment?\n 1 = Simple Correlation\n 2 = Haywood\n");
AF_selection = input("Autofocus?\n 1 = Yuan \n 2 = Haywood\n");
%% 2 Simulation Parameters
C = 3e8;                                            % Speed of light
F0 = 9.5e9;                                         % Initial centre frequency
DeltaF = 4e6;                                       % Frequency step size
N = 256;                                            % Number of pulses in burst
M = 32;                                             % Number of Bursts 
PRF = 2e3;                                          % Pulse Repetition Frequency
SNR_dB = 40;                                        % Assume signal power = 1
CentreFrequencyVector = (0:1:(N-1))*DeltaF + F0;
R0 = 10.106e3;                                      % Distance to centre of target in km 
Pn = 10^(-(SNR_dB-10*log10(N))/10);                 % Noise Power

Range_Resolution = C/(2*N*DeltaF);
Unambiguous_range = (N-1)*Range_Resolution;
PRI = 1/PRF;
BurstRepetionFrequency = PRF/M; 
BurstRepetionInterval = 1/BurstRepetionFrequency;

TimePerBurstVector = (0:1:(M-1))*BurstRepetionInterval;

disp(' ');
disp(['Range resolution = ' num2str(roundn(Range_Resolution,-3)) ' m']);
disp(['Unambiguous Range = ' num2str(roundn(Unambiguous_range,0)) ' m']);
disp(' ');

%% 3 Define Scatterers and Target motion

TgtVelocity_ms = target_velocity;                   % Translational motion
RotRate_deg_s = rotation_angle;                     % Rotation angle

% Scatterer local co-ordinates
% Each row has two elements: [x-cord y-cord]  
                        % x y amplitude
Scatterer_axy_local = [
    -10 0 1; -9 0 1; -8 0 1; -7 0 1.5; -6 0 2; -5 0 2.5; -4 0 2.5; 
    -3 0 3;-3 1 3;-3 2 2;-3 3 1;
    -2 0 5; -1 0 6; 0 0 20; 
    1 0 6; 2 0 5; 
    3 0 3; 3 1 3;3 2 2;3 3 1;
    4 0 2.5; 5 0 2.5; 6 0 2; 7 0 1.5; 8 0 1; 9 0 1; 10 0 1];


% Calculate the distance of each point from the center of the object (0,0)
distances = sqrt(sum(Scatterer_axy_local(:, 1:2).^2, 2));

% Calculate the amplitude values based on a Gaussian-like distribution
% The centermost scatterer will have the highest amplitude
max_amplitude = 20; % Set the maximum amplitude
std = 5; % Adjust this parameter to control the distribution width
amplitudes = max_amplitude * exp(-(distances.^2) / (2 * std^2));

% Replace the amplitude values in the scatterer data
Scatterer_axy_local(:, 3) = amplitudes;

figure; scatter(Scatterer_axy_local(:,1),Scatterer_axy_local(:,2))
xlabel('x-coordinate (m)');
ylabel('y-coordinate');
title('Plot of the target scatterers');
%% 4 Simulation
NumScatterers = size(Scatterer_axy_local,1);               % Obtain number of scatterers

for BurstCounter = 1:M
    
    RotationAngleDeg(BurstCounter) = RotRate_deg_s*(BurstCounter-1)*BurstRepetionInterval;
    
     for ScattererCounter = 1:NumScatterers
         % Rotate scatterers one at a time, for every burst
         SingleScatterer_xy_local_Original = Scatterer_axy_local(ScattererCounter, 1:2).'; % get the xy-coords
         RotationMatrix = [ cosd(RotationAngleDeg(BurstCounter)) -sind(RotationAngleDeg(BurstCounter)); sind(RotationAngleDeg(BurstCounter)) cosd(RotationAngleDeg(BurstCounter))];
         SingleScatterer_xy_local_rotated = RotationMatrix*SingleScatterer_xy_local_Original;
         Scatterer_xy_local_rotated(ScattererCounter, :) = SingleScatterer_xy_local_rotated';
     end
    
    % Define where the scatterers are in terms of the global axis (axis of
    % radar) - R0 is distance from radar to center of target
    Scatterer_xy_global(:,1) = Scatterer_xy_local_rotated(:,1) + R0 + TgtVelocity_ms*BurstRepetionInterval*(BurstCounter-1);
    Scatterer_xy_global(:,2) = Scatterer_xy_local_rotated(:,2);
    
    % calculate the distance from the radar to the scatterer
    Rk = sqrt(Scatterer_xy_global(:,1).^2 + Scatterer_xy_global(:,2).^2);   % Slant range to each Scatterer
    Rx = zeros(1,N);
    
    % For all n:
    % Rx(n) = An * exp(-1i*2*pi*CentreFrequencyVector*2*Rk(n)/C) 
    %       = An * exp(-1i * 2 * pi * fc * tau) ; tau = 2Rk(n)/C
    for CountScatterer = 1:  NumScatterers
        Rx = Rx + Scatterer_axy_local(CountScatterer,3)*exp(-1i*4*pi*CentreFrequencyVector*Rk(CountScatterer)/C);
    end
    
    %% 4.1 Add noise
    a = randn(1,N);
    b = randn(1,N);
    
    NoiseVoltageVector = sqrt(Pn)*(a+CountScatterer*b)*1/sqrt(2);  %
    RxNoise(BurstCounter, :) = Rx + NoiseVoltageVector;                                    % Receive Signal plus noise
    
end

%% 5 Get Unfocused HRR Profiles
HRR_profiles = (ifft(RxNoise,[], 2));
RangeAxis = (0:1:(N-1))*C/(2*N*DeltaF); 

% Plot HRR profiles
figure;
linear_dB = Normalise_limitDynamicRange_ISAR_dB(HRR_profiles,35);
imagesc(RangeAxis, 1:M, linear_dB);
xlabel('Range (m)');
ylabel('Profile Number');
title('Unaligned HRR Profiles');
colormap('jet');
colorbar;

% Plot unfocused ISAR image
WindowMatrix = repmat(hamming(M),1, N);
ISAR_linear = fftshift(fft(HRR_profiles.*WindowMatrix, [], 1),1);
ISAR_linear_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_linear,40);
FrequencyAxis_Hz = (-M/2:1:(M/2-1))*BurstRepetionFrequency/M;

figure;
imagesc(RangeAxis, FrequencyAxis_Hz, ISAR_linear_dB);
xlabel('Range (m)');
ylabel('Doppler frquency (Hz)');
title('Unfocused ISAR image');
colorbar;
colormap('jet');
axis xy;
%% 6 Range Alignment of Profiles 
% Range Align the HRR profiles using user selected algorithm
ref_profile_number =1;
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

% Plot HRR profiles
figure;
linear_dB = Normalise_limitDynamicRange_ISAR_dB(RA_HRR_profiles,35);
imagesc(RangeAxis, 1:M, linear_dB);
xlabel('Range (m)');
ylabel('Profile Number');
title('Haywood Range-aligned HRR Profiles');
colormap('jet');
colorbar;

% Plot ISAR image
WindowMatrix = repmat(hamming(M),1, N);
ISAR_linear = fftshift(fft(RA_HRR_profiles.*WindowMatrix, [], 1),1);
ISAR_linear_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_linear,35);
FrequencyAxis_Hz = (-M/2:1:(M/2-1))*BurstRepetionFrequency/M;

figure;
imagesc(RangeAxis, FrequencyAxis_Hz, ISAR_linear_dB);
xlabel('Range (m)');
ylabel('Doppler frquency (Hz)');
title('Haywood Range-aligned ISAR image');
colorbar;
colormap('jet');
axis xy;

%% Autofocus of Profiles
% Autofocus the HRR profiles using selected method
if(AF_selection == 1)
    [AF_RA_HRR_profiles] = YuanAF(RA_HRR_profiles);
elseif (AF_selection == 2)
    [AF_RA_HRR_profiles] = HaywoodAF(RA_HRR_profiles);

end

% Plot ISAR image
WindowMatrix = repmat(hamming(M),1, N);
ISAR_linear = fftshift(fft(AF_RA_HRR_profiles.*WindowMatrix, [], 1),1);
ISAR_linear_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_linear,35);
FrequencyAxis_Hz = (-M/2:1:(M/2-1))*BurstRepetionFrequency/M;

figure;
imagesc(RangeAxis, FrequencyAxis_Hz, ISAR_linear_dB);
xlabel('Range (m)');
ylabel('Doppler frquency (Hz)');
title('Focused ISAR image');
colorbar;
colormap('jet');
axis xy;
