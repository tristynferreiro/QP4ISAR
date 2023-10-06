clear all;
close all;
clc;


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

Scatterer_xy_local = [-10 0; 0 0; 7 0 ];                        % Scatterer local co-ordinates
                                                          % Each row has two elements: [x-cord y-cord]  
NumScatterers = size(Scatterer_xy_local,1);               % Obtain number of scatterers
Scatterer_xy_global(:,1) = Scatterer_xy_local(:,1) + R0;
Scatterer_xy_global(:,2) = Scatterer_xy_local(:,2);

Rk = sqrt(Scatterer_xy_global(:,1).^2 + Scatterer_xy_global(:,2).^2);   % Slant range to each Scatterer
Rx = zeros(1,N);

for CountScatterer = 1:  NumScatterers
 Rx = Rx + exp(-1i*4*pi*CentreFrequencyVector*Rk(CountScatterer)/C);
end

a = randn(1,N);
b = randn(1,N);

NoiseVoltageVector = sqrt(Pn)*(a+i*b)*1/sqrt(2);  % 
RxN = Rx + NoiseVoltageVector;                                    % Receive Signal plus noise

HRR_Profile = (ifft(RxN));
RangeAxis = (0:1:(N-1))*Range_Resolution;

figure;
plot(RangeAxis, 20*log(abs(HRR_Profile)))
grid on;

x = 1;


