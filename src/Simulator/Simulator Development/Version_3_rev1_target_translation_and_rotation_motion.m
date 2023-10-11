%% Version 3 - include target translation and rotational motion 
% Revision: allow configurable scatterer amplitudes.
clear all; close all; clc;

%% Simulation Parameters
C = 3e8;                                            % Speed of light
F0 = 9.5e9;                                         % Initial centre frequency
DeltaF = 4e6;                                       % Frequency step size
N = 256;                                            % Number of pulses in burst
M = 32;                                             % Number of Bursts 
PRF = 2e3;                                          % Pulse Repetition Frequency
SNR_dB = 40;                                        % Assume signal power = 1
CentreFrequencyVector = (0:1:(N-1))*DeltaF + F0;
R0 = 10.106e3;                                           % Distance to centre of target in km 
Pn = 10^(-(SNR_dB-10*log10(N))/10);                % Noise Power

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

TgtVelocity_ms = -1; 
RotRate_deg_s = 6;
% Scatterer local co-ordinates
% Each row has two elements: [x-cord y-cord]  
                    % x y amplitude
Scatterer_axy_local = [-10 0 1; -9 0 1; -8 0 1; -7 0 1.5; -6 0 2; -5 0 2.5; -4 0 2.5; -3 0 3; -2 0 5; -1 0 6; 0 0 20; 
                         1 0 6; 2 0 5; 3 0 3; 4 0 2.5; 5 0 2.5; 6 0 2; 7 0 1.5; 8 0 1; 9 0 1; 10 0 1; 
-3 1 3; 3 1 3; -3 2 2; 3 2 2; -3 3 1; 3 3 1];

figure; scatter(Scatterer_axy_local(:,1),Scatterer_axy_local(:,2))

%% Target scatterers
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
    
    %% Add noise
    a = randn(1,N);
    b = randn(1,N);
    
    NoiseVoltageVector = sqrt(Pn)*(a+CountScatterer*b)*1/sqrt(2);  %
    RxNoise(BurstCounter, :) = Rx + NoiseVoltageVector;                                    % Receive Signal plus noise
    
end

%% Get HRR Profiles
HRR_Profile = (ifft(RxNoise,[], 2));
RangeAxis = (0:1:(N-1))*C/(2*N*DeltaF); 

% Plot HRR profiles
figure;
linear_dB = Normalise_limitDynamicRange_ISAR_dB(HRR_Profile,35);
imagesc(RangeAxis, 1:M, linear_dB);
xlabel('Range (m)');
ylabel('Profile Number');
title('Unaligned HRR Profiles');
colormap('jet');
colorbar;

% Plot unfocused ISAR image
WindowMatrix = repmat(hamming(M),1, N);
ISAR_linear = fftshift(fft(HRR_Profile.*WindowMatrix, [], 1),1);
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