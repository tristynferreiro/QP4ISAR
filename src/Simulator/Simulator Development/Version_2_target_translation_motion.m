
% Version 2 - include target translation motion 

clear all;
close all;
clc;

C = 3e8;                                            % Speed of light
F0 = 9.5e9;                                         % Initial centre frequency
DeltaF = 4e6;                                       % Frequency step size
N = 256;                                            % Number of pulses in burst
M = 160;                                             % Number of Bursts 
PRF = 6e3;                                          % Pulse Repetition Frequency
SNR_dB = 40;                                        % Assume signal power = 1
CentreFrequencyVector = (0:1:(N-1))*DeltaF + F0;
R0 = 10e3;                                           % Distance to centre of target in km 
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

TgtVelocity_ms = 10; 
Scatterer_axy_local = [-10 0 30; -7 0 20; -3 0 40; 0 0 60; 3 0 40; 7 0 20; 10 0 30; 
    0 1 10; 0 2 10; 0 3 10];              
                                                          % Scatterer local co-ordinates
                                                          % Each row has two elements: [x-cord y-cord]  
NumScatterers = size(Scatterer_axy_local,1);               % Obtain number of scatterers

for BurstCounter = 1:M
    
    Scatterer_xy_global(:,1) = Scatterer_axy_local(:,1) + R0 + TgtVelocity_ms*BurstRepetionInterval*(BurstCounter);
    Scatterer_xy_global(:,2) = Scatterer_axy_local(:,2);
    
    Rk = sqrt(Scatterer_xy_global(:,1).^2 + Scatterer_xy_global(:,2).^2);   % Slant range to each Scatterer
    Rx = zeros(1,N);
    
    for CountScatterer = 1:  NumScatterers
        Rx = Rx + Scatterer_axy_local(CountScatterer,3)*exp(-1i*4*pi*CentreFrequencyVector*Rk(CountScatterer)/C);
    end
    
    a = randn(1,N);
    b = randn(1,N);
    
    PnVector = sqrt(Pn)*(a+1i*b)*1/sqrt(2);  %
    RxNoise(BurstCounter, :) = Rx + PnVector;                                    % Receive Signal plus noise
    
end

HRR_Profile = (ifft(RxNoise,[], 2));
RangeAxis = (0:1:(N-1))*C/(2*N*DeltaF);

figure;
HRR_Profile_dB = Normalise_limitDynamicRange_ISAR_dB(HRR_Profile,SNR_dB);
imagesc(RangeAxis, 1:M, HRR_Profile_dB);
xlabel('Range (m)');
ylabel('Profile Number');
title('HRR Profiles');
%colormap('jet');
colorbar;
grid on;

% Plot ISAR image
WindowMatrix = repmat(hamming(M),1, N);
ISAR_linear = fftshift(fft(HRR_Profile.*WindowMatrix, [], 1),1);
ISAR_linear_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_linear,SNR_dB);
FrequencyAxis_Hz = (-M/2:1:(M/2-1))*BurstRepetionFrequency/M;

figure;
imagesc(RangeAxis, FrequencyAxis_Hz, ISAR_linear_dB);
xlabel('Range (m)');
ylabel('Doppler frquency (Hz)');
title('Unfocused ISAR image');
colorbar;
colormap('jet');
axis xy;

%% Range Alignment of Profiles 
% Range Align the HRR profiles using correlation method
ref_profile_number =1;
[corr_RA_HRR_profiles] = correlationRA(HRR_Profile,ref_profile_number);
%[haywood_RA_HRR_profiles] = HaywoodRA(HRR_profiles,ref_profile_number);

% Plot HRR profiles
figure;
HRR_Profile_dB = Normalise_limitDynamicRange_ISAR_dB(corr_RA_HRR_profiles,SNR_dB);
imagesc(RangeAxis, 1:M, HRR_Profile_dB);
xlabel('Range (m)');
ylabel('Profile Number');
title('HRR Profiles');
%colormap('jet');
colorbar;
grid on;

% Plot ISAR image
WindowMatrix = repmat(hamming(M),1, N);
ISAR_linear = fftshift(fft(corr_RA_HRR_profiles.*WindowMatrix, [], 1),1);
ISAR_linear_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_linear,SNR_dB);
FrequencyAxis_Hz = (-M/2:1:(M/2-1))*BurstRepetionFrequency/M;

figure;
imagesc(RangeAxis, FrequencyAxis_Hz, ISAR_linear_dB);
xlabel('Range (m)');
ylabel('Doppler frquency (Hz)');
title('ISAR image');
colorbar;
colormap('jet');
axis xy;

%% Autofocus of Profiles
%Apply Haywood autofocus to the RA HRR profiles using
AF_corrRA_HRR_profiles = HaywoodAF(corr_RA_HRR_profiles);

% Plot ISAR image
WindowMatrix = repmat(hamming(M),1, N);
ISAR_linear = fftshift(fft(AF_corrRA_HRR_profiles.*WindowMatrix, [], 1),1);
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



