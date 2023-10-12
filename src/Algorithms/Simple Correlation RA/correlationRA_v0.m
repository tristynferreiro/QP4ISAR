function [RA_HRRP,shifts] = correlationRA_v0(HRRP_all, ref_profile_number)
% Implements a simple correlation range-alignment algorithm. 
    % This approach aligns all profiles based on the a selected reference 
    % profile. We align the profiles based on the peak values in each 
    % profile and use correlation to discover the location of the peaks. 
    % This location corresponds to the number of range bins the profile 
    % needs to be shifted in order to align with the reference profile.
    
    % Principles of Modern Radar was used as a resource in
    % implementing this algorithm
    %% Parameters
    RA_HRRP = HRRP_all;
    shifts = zeros(size(HRRP_all,1),1);
    
    %% Step 1: Autocorrelation, using the reference profile
    ref_HRRP= HRRP_all(ref_profile_number,:);
    correlationOut = xcorr(abs(ref_HRRP), abs(ref_HRRP));
    PeakIndexRef = find(correlationOut==max(correlationOut));
    
    %% Step 2: Range Alignment using ref profile
    % Now begin alignment with all other profiles:
    [rows] = size(HRRP_all);
    
    for indx = 1:rows % profile 1 is the reference
        rangeProfile = HRRP_all(indx,:); % Get range profile
        % Calculate correlation
        correlationOut = xcorr(abs(ref_HRRP), abs(rangeProfile));
        % Find index of peak correlation value
        PeakIndex = find(correlationOut==max(correlationOut));
        % Calculate shift between ref peak and this range profile peak
        shift = PeakIndex - PeakIndexRef;
        shifts(indx)=shift;

        % Shift range profile to align with ref
        shiftedRangeProfile = circshift(rangeProfile,shift);
        % Add shifted profile to array
        RA_HRRP(indx,:) = shiftedRangeProfile;
    end
end