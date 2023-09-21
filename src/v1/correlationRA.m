function [RA_HRRP,shifts] = correlationRA(HRRP, ref_HRRP_num)
% Implements a simple correlation range-alignment algorithm. 
    % This approach aligns all profiles based on the a selected reference 
    % profile. We align the profiles based on the peak values in each 
    % profile and use correlation to discover the location of the peaks. 
    % This location corresponds to the number of range bins the profile 
    % needs to be shifted in order to align with the reference profile.
    
    %% Step 1 : get reference profile
    ref_HRRP= HRRP(ref_HRRP_num,:);
   
    %% Step 2: Compute the shift values between all profiles and the ref
    % compute correlation values
    correlation = xcorr2(abs(HRRP),abs(ref_HRRP));

     % Find index of peak correlation values
    [~,peak_index] = max(correlation,[],2);

    % Calculate shifs between ref peak and all profile peaks
    shifts = peak_index(1) - peak_index;
    
    %% Step 3: Range Alignment using ref profile
    %RA_HRRP = circshift(HRRP,shifts');
    RA_HRRP = HRRP;
    rows = size(HRRP,1);
    for indx = 1:rows
        % Shift range profile to align with ref and add to array
        RA_HRRP(indx,:) = circshift(HRRP(indx,:),shifts(indx));
    end

end