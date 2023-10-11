
% Version 2 - include target translation motion 

clear all;
close all;
clc;

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

TgtVelocity_ms = 0; 
RotRate_deg_s = 6;
%Scatterer_xy_local = [-10 0; 0 0; 7 0 ];                  % Scatterer local co-ordinates
                                                          % Each row has two elements: [x-cord y-cord]  

Scatterer_xy_local = [-10 0; -7 0; -3 0; 0 0; 3 0; 7 0; 10 0; -3 1; 3 1; -3 2; 3 2; -3 3; 3 3 ];
                                                          
NumScatterers = size(Scatterer_xy_local,1);               % Obtain number of scatterers

for BurstCounter = 1:M
    
    RotationAngleDeg(BurstCounter) = RotRate_deg_s*(BurstCounter-1)*BurstRepetionInterval;
    
     for ScattererCounter = 1:NumScatterers
         % Rotate scatterers one at a time, for every burst
         
         SingleScatterer_xy_local_Original = Scatterer_xy_local(ScattererCounter, :).';
         RotationMatrix = [ cosd(RotationAngleDeg(BurstCounter)) -sind(RotationAngleDeg(BurstCounter)); sind(RotationAngleDeg(BurstCounter)) cosd(RotationAngleDeg(BurstCounter))];
         SingleScatterer_xy_local_rotated = RotationMatrix*SingleScatterer_xy_local_Original;
         Scatterer_xy_local_rotated(ScattererCounter, :) = SingleScatterer_xy_local_rotated';
     end
    
    Scatterer_xy_global(:,1) = Scatterer_xy_local_rotated(:,1) + R0 + TgtVelocity_ms*BurstRepetionInterval*(BurstCounter-1);
    Scatterer_xy_global(:,2) = Scatterer_xy_local_rotated(:,2);
    
    Rk = sqrt(Scatterer_xy_global(:,1).^2 + Scatterer_xy_global(:,2).^2);   % Slant range to each Scatterer
    Rx = zeros(1,N);
    
    for i = 1:  NumScatterers
        Rx = Rx + exp(-1i*4*pi*CentreFrequencyVector*Rk(i)/C);
    end
    
    a = randn(1,N);
    b = randn(1,N);
    
    PnVector = sqrt(Pn)*(a+i*b)*1/sqrt(2);  %
    RxN(BurstCounter, :) = Rx + PnVector;                                    % Receive Signal plus noise
    
end

HRR_Profile = (ifft(RxN,[], 2));
RangeAxis = (0:1:(N-1))*C/(2*N*DeltaF);

% Plot HRR profiles

figure;
HRR_Profile_dB = Normalise_limitDynamicRange_ISAR_dB(HRR_Profile,40);
imagesc(RangeAxis, 1:M, HRR_Profile_dB);
xlabel('Range (m)');
ylabel('Profile Number');
title('HRR Profiles');
%colormap('jet');
colorbar;


% Plot ISAR image
WindowMatrix = repmat(hamming(M),1, N);;
ISAR_linear = fftshift(fft(HRR_Profile.*WindowMatrix, [], 1),1);
ISAR_linear_dB = Normalise_limitDynamicRange_ISAR_dB(ISAR_linear,40);
FrequencyAxis_Hz = (-M/2:1:(M/2-1))*BurstRepetionFrequency/M;

figure;
imagesc(RangeAxis, FrequencyAxis_Hz, ISAR_linear_dB);
xlabel('Range (m)');
ylabel('Doppler frquency (Hz)');
title('ISAR image');
colorbar;
colormap('jet');
axis xy;

x =1;

