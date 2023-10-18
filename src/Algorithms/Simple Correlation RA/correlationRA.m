function [RA_HRRP] = correlationRA(HRRP_all, ref_HRRP_num)
% Implements a simple correlation range-alignment algorithm. 
    % This approach aligns all profiles based on the a selected reference 
    % profile. We align the profiles based on the peak values in each 
    % profile and use correlation to discover the location of the peaks. 
    % This location corresponds to the number of range bins the profile 
    % needs to be shifted in order to align with the reference profile.
    %
    % Principles of Modern Radar was used as a resource in
    % implementing this algorithm
    % 
    % Version: v1.1
    %
    % Revisions:
    %   - Further optimisation changes using built in MATLAB 
    %     functions, reducing calcualtions within for loops and general 
    %     code cleaning
    
    %% Step 1 : get reference profile
    ref_HRRP= HRRP_all(ref_HRRP_num,:);
   
    %% Step 2: Compute the shift values between all profiles and the ref
    % Compute correlation values
    correlation = xcorr2(abs(HRRP_all),abs(ref_HRRP));

     % Find index of peak correlation values
    [~,peak_index] = max(correlation,[],2);

    % Calculate shifs between ref peak and all profile peaks
    shifts = peak_index(ref_HRRP_num) - peak_index;
    % figure; plot(peak_index(1)-peak_index-shifts) % check shifts are
    % correct
    
    % Plot stair case Function
    % figure; plot(1:size(shifts,1),shifts)
    % xlabel('Profile Number');
    % ylabel('Number of bin shifts')
    % title('Bin shifts per Range Profile');

    %% Step 3: Range Alignment using ref profile
    RA_HRRP = HRRP_all; % pre-defined for efficiency
    num_profiles = size(HRRP_all,1);
    
    for indx = 1:num_profiles
        % Shift range profile to align with ref and add to array
        RA_HRRP(indx,:) = circshift(HRRP_all(indx,:),shifts(indx));
    end

end