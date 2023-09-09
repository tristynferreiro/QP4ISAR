%% Implementation of the Haywood RA algorithm
% Initial explanation for understanding:
RangeProfiles = [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 
    0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0];
figure; imagesc(abs(RangeProfiles))
colorbar
%% Simple Test Case
HRRP_ref = [0 0 1 0 0 0 0 0];
d = 0.5; 
n=8;
m = 0:1:n-1;

phi = exp(-1i*(2*pi*d*m)/n);
shifted_test = ifft(phi.*fft(HRRP_ref));
disp(abs(shifted_test))

figure; imagesc(abs(shifted_test))
colorbar
%% Real implementation
% This uses the Haywood Range-Alignment algorithm (no autofocus). Haywood's
% algorithm begins with the simple correlation
ref_profile_number =4;

% Step 1: Get shift values
[corr_HRR_profiles,shifts] = correlationRA(RangeProfiles, ref_profile_number);

% Plot Step Function
figure; plot(1:size(shifts,1),shifts)
xlabel('Profile Number');
ylabel('Number of bin shifts')
title('Bin shifts per Range Profile');
hold on

% Step 2: smooth the shift values so we polyfit to the lowest degree
% polynomial
numProfiles = 1:size(shifts);
coefficients = polyfit(numProfiles,shifts,1);
s = polyval(coefficients,numProfiles); % smoothed values

% plot smoothed values
plot(numProfiles,s,'-')
hold off

% Step 3: range-align profiles
% shifted=RangeProfiles;

% for i = 1: size(s,2)
%     phi = exp(-1i*(2*pi*s(i)*m)/n);
%     %disp(RangeProfiles(i,:))
%     shifted(i,:) = ifft(phi.*fft(RangeProfiles(i,:)));
%     %disp(abs(shifted))
% end

% figure; imagesc(abs(range))
% colorbar
% figure; imagesc(abs(shifted))
% colorbar
%% 

