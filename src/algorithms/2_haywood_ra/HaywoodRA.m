function [RA_HRRP] = HaywoodRA(HRRP_all, ref_HRRP_num)
% Implements the Haywood range-alignment algorithm
    % This approach aligns all profiles based on the a selected reference 
    % profile. We align the profiles based on the peak values in each 
    % profile and use correlation to discover the location of the peaks. 
    % This location corresponds to the number of range bins the profile 
    % needs to be shifted by. These shift values are then linearised and 
    % used to calculate a phase adjustment factor for each profile. The 
    % phase adjustment factor is then applied to each profile in order to 
    % align with the reference profile.

    % A. Zyweck, PhD Thesis Appendix was used as a resource in
    % implementing this Haywood algorithm

    % Revision 2: Implements cross-correlation within the script, rather
    % than calling another script

    %% Step 1 get the correlation values (the bin shifts): - Eq A.3 of Zyweck's appendix
    % Get reference profile
    ref_HRRP= HRRP_all(ref_HRRP_num,:);
   
    % compute correlation values
    correlation = xcorr2(abs(HRRP_all),abs(ref_HRRP));

    % Find index of peak correlation values
    [~,peak_index] = max(correlation,[],2);

    % Calculate shifs between ref peak and all profile peaks
    shifts = peak_index(ref_HRRP_num) - peak_index;
    % figure; plot(peak_index(1)-peak_index-shifts) % check shifts are
    % correct
    
    % % Plot stair case Function
    % figure; plot(1:size(shifts,1),shifts)
    % xlabel('Profile Number');
    % ylabel('Number of bin shifts')
    % title('Bin shifts per Range Profile');
    % hold on

    %% Step 2: smooth the shift values so we polyfit to the lowest degree
    % general parameters
    [num_profiles,num_range_bins]=size(HRRP_all);
    profiles = 1:num_profiles;

    % linearise the shift values
    coefficients = polyfit(profiles,shifts,1);
    smoothed_coefficients = polyval(coefficients,profiles); % smoothed values

    % Plot smoothed values
    % plot(profiles,smoothed_coefficients,'-')
    % hold off
    % matlab2tikz() % Only uncomment to Save the figure as LaTeX compatible plot

    %% Step 3: Perform Haywood's Range alignment

    % calculate phase correction factor - Eq A.6 of Zyweck's appendix
    m = 0:1:num_range_bins-1;
    phase = exp(-1i*(2*pi*smoothed_coefficients'.*m)/num_range_bins);

    % Apply correction angle to profiles - Eq A.5 of Zyweck's appendix
    RA_HRRP = ifft(phase.*fft(HRRP_all,num_range_bins,2),num_range_bins,2);
end