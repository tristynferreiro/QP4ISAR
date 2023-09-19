function [RA_HRR_profiles] = HaywoodRA(HRR_profiles, ref_profile_number)
% Implements the Haywood range-alignment algorithm
    % This approach uses the correlation RA algorithm but incorporates
    % phase adjustment to achieve better range-alignment.


%% Step 1 get the correlation values (the bin shifts):
[~,shifts] = correlationRA(HRR_profiles,ref_profile_number);

% Plot stair case Function
% figure; plot(1:size(shifts,1),shifts)
% xlabel('Profile Number');
% ylabel('Number of bin shifts')
% title('Bin shifts per Range Profile');
% hold on

%% Step 2: smooth the shift values so we polyfit to the lowest degree
% polynomial
[num_profiles,num_range_bins]=size(HRR_profiles);
profiles = 1:num_profiles;
coefficients = polyfit(profiles,shifts,1);
smoothed_coefficients = polyval(coefficients,profiles); % smoothed values

% plot smoothed values
% plot(numProfiles,smoothed_coefficients,'-')
% hold off

%% Step 3: Perform Haywood's Range alignment
m = 0:1:num_range_bins-1;
% calculate phase correction angle
phase = exp(-1i*(2*pi*smoothed_coefficients'.*m)/num_profiles);
% Apply correction angle to profiles
RA_HRR_profiles = ifft(phase.*fft(HRR_profiles,num_range_bins,2),num_range_bins,2);
end