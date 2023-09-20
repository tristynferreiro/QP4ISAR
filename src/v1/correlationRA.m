function [RA_HRRP,shifts] = correlationRA(HRRP, ref_HRRP_num)
% Implements a simple correlation range-alignment algorithm. 
    % This approach aligns all profiles based on the a selected reference 
    % profile. We align the profiles based on the peak values in each 
    % profile and use correlation to discover the location of the peaks. 
    % This location corresponds to the number of range bins the profile 
    % needs to be shifted in order to align with the reference profile.

    RA_HRRP = HRRP;
    shifts = zeros(size(HRRP,1),1);

    %% Step 1: Autocorrelation, using the reference profile
    ref_HRRP= HRRP(ref_HRRP_num,:);
    corr_out = xcorr(abs(ref_HRRP));
    [~,peak_index_ref] = max(corr_out);
    
    %% Step 2: Range Alignment using ref profile
    rows = size(HRRP);

    for indx = 1:rows
        % Get range profile
        range_profile = HRRP(indx,:); 
        
        % Calculate correlation
        corr_out = xcorr(abs(ref_HRRP), abs(range_profile));
        
        % Find index of peak correlation value
        [~,peak_index] = max(corr_out);
      
        % Calculate shift between ref peak and this range profile peak
        shift = peak_index - peak_index_ref;
        shifts(indx)=shift;

        % Shift range profile to align with ref and add to array
        RA_HRRP(indx,:) = circshift(range_profile,shift);
    end
end