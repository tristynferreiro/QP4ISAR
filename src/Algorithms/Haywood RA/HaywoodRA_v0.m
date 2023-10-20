function [RA_HRRP] = HaywoodRA_v0(HRRP_all, ref_profile_number)
% Implements the Haywood range-alignment algorithm
    % This approach aligns all profiles based on the a selected reference 
    % profile. We align the profiles based on the peak values in each 
    % profile and use correlation to discover the location of the peaks. 
    % This location corresponds to the number of range bins the profile 
    % needs to be shifted by. These shift values are then linearised and 
    % used to calculate a phase adjustment factor for each profile. The 
    % phase adjustment factor is then applied to each profile in order to 
    % align with the reference profile.
    %
    % A. Zyweck, PhD Thesis, Appendix was used as a resource in
    % implementing this Haywood algorithm
    %
    % Version: v0

%% Step 1 get the correlation values (the bin shifts):
[RA_HRRP,shifts] = correlationRA_v0(HRRP_all,ref_profile_number);

% Plot stair case Function
% figure; plot(1:size(shifts,1),shifts)
% xlabel('Profile Number');
% ylabel('Number of bin shifts')
% title('Bin shifts per Range Profile');
% hold on

%% Step 2: smooth the shift values so we polyfit to the lowest degree
% polynomial
profile_nums = 1:size(shifts);
coefficients = polyfit(profile_nums,shifts,1);
s = polyval(coefficients,profile_nums); % smoothed values

% plot smoothed values
% plot(numProfiles,s,'-')
% hold off

%% Step 3: Perform Haywood's Range alignment
[profiles,bins]=size(HRRP_all);

m = 0:1:bins-1;

for i = 1: profiles
    phi = exp(-1i*(2*pi*s(i)*m)/bins);  % Calculate phase shift
    %disp(RangeProfiles(i,:))
    RA_HRRP(i,:) = ifft(phi.*fft(HRRP_all(i,:))); % shift profiles
    %disp(abs(shifted))
end

end