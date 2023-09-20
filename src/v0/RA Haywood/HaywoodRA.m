function [RA_HRR_profiles] = HaywoodRA(HRR_profiles, ref_profile_number)
% Implements the Haywood range-alignment algorithm
    % This approach uses the correlation RA algorithm but incorporates
    % phase adjustment to achieve better range-alignment.

%% Step 1 get the correlation values (the bin shifts):
[RA_HRR_profiles,shifts] = correlationRA(HRR_profiles,ref_profile_number);

% Plot stair case Function
% figure; plot(1:size(shifts,1),shifts)
% xlabel('Profile Number');
% ylabel('Number of bin shifts')
% title('Bin shifts per Range Profile');
% hold on

%% Step 2: smooth the shift values so we polyfit to the lowest degree
% polynomial
numProfiles = 1:size(shifts);
coefficients = polyfit(numProfiles,shifts,1); % fit to straight line
s = polyval(coefficients,numProfiles); % smoothed values

% plot smoothed values
% plot(numProfiles,s,'-')
% hold off

%% Step 3: Perform Haywood's Range alignment
[profiles,bins]=size(HRR_profiles);

m = 0:1:bins-1;

for i = 1: profiles
    phi = exp(-1i*(2*pi*s(i)*m)/profiles);
    %disp(RangeProfiles(i,:))
    RA_HRR_profiles(i,:) = ifft(phi.*fft(HRR_profiles(i,:)));
    %disp(abs(shifted))
end

end