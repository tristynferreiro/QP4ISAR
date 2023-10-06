clear all; close all; clc;

% This is an updated version of the original code with the following 
% changes:
%   1. The amplitude of scatterers can be configured (previously only the
%      x-y coords could be set)
%   2. Additional comments
%% Simulation Parameters
C = 3e8;                                            % Speed of light
F0 = 9.5e9;                                         % Initial centre frequency
DeltaF = 4e6;                                       % Frequency step size
N = 256;                                            % Number of pulses in burst
SNR_dB = 80;                                        % Assume signal power = 1
CentreFrequencyVector = (0:1:(N-1))*DeltaF + F0;
R0 = 10e3;                                           % Distance to centre of target in km 
Pn = 10^(-(SNR_dB/2-10*log10(N))/10);                % Noise Power

Range_Resolution = C/(2*N*DeltaF);
Unambiguous_range = (N-1)*Range_Resolution;

disp(' ');
disp(['Range resolution = ' num2str(roundn(Range_Resolution,-3)) ' m']);
disp(['Unambiguous Range = ' num2str(roundn(Unambiguous_range,0)) ' m']);
disp(' ');

%% Target scatterers
                     % amplitude x y
Scatterer_axy_local = [1 -10 0; 2 0 0; 1 7 0 ];               % Scatterer amplitude and local co-ordinates
                                                       % Each row (scatterer) has two elements: [x-cord y-cord]  
NumScatterers = size(Scatterer_axy_local,1);            % Obtain number of scatterers

% Define where the scatterers are in terms of the global axis (axis of
% radar) - R0 is distance from radar to center of target
Scatterer_xy_global(:,1) = Scatterer_axy_local(:,2) + R0;
Scatterer_xy_global(:,2) = Scatterer_axy_local(:,3);

% calculate the distance from the radar to the scatterer
Rk = sqrt(Scatterer_xy_global(:,1).^2 + Scatterer_xy_global(:,2).^2);   % Slant range to each Scatterer
Rx = zeros(1,N);

% For all n:
% Rx(n) = An * exp(-1i*2*pi*CentreFrequencyVector*2*Rk(n)/C) 
%       = An * exp(-1i * 2 * pi * fc * tau) ; tau = 2Rk(n)/C
for CountScatterer = 1:  NumScatterers
 Rx = Rx + Scatterer_axy_local(CountScatterer,1)*exp(-1i*4*pi*CentreFrequencyVector*Rk(CountScatterer)/C);
end

%% Add noise
a = randn(1,N);
b = randn(1,N);

NoiseVoltageVector = sqrt(Pn)*(a+1i*b)*1/sqrt(2);  % 
RxNoise = Rx + NoiseVoltageVector;                                    % Receive Signal plus noise

%% Get HRR Profiles
HRR_Profile = (ifft(RxNoise));
RangeAxis = (0:1:(N-1))*Range_Resolution;

figure;
plot(RangeAxis, 20*log(abs(HRR_Profile)))
xlabel('Range (m)');
ylabel('Voltage (dB)');
title('HRR Profiles');
grid on;



