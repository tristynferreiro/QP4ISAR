function [RA_HRR_profiles,shifts] = correlationRA(HRR_profiles, ref_profile_number)
% Implements a simple correlation range-alignment algorithm. 
    % This approach aligns all profiles based on the a selected reference 
    % profile. We align the profiles based on the peak values in each 
    % profile and use correlation to discover the location of the peaks. 
    % This location corresponds to the number of range bins the profile 
    % needs to be shifted in order to align with the reference profile.

    RA_HRR_profiles = HRR_profiles;
    shifts = zeros(size(HRR_profiles,1),1);

    % Autocorrelation, choosing the 1st profile to be the reference profile:
    refHRRProfile= HRR_profiles(ref_profile_number,:);
    CorrelationOut = xcorr(abs(refHRRProfile), abs(refHRRProfile));
    PeakIndexRef = find(CorrelationOut==max(CorrelationOut));
    %RA_HRR_profiles(ref_profile_number,:) = refHRRProfile;
    
    % Now begin alignment with all other profiles:
    [rows] = size(HRR_profiles);
    for indx = 1:rows % profile 1 is out reference
        rangeProfile = HRR_profiles(indx,:); % Get range profile
        % Calculate correlation
        CorrelationOut = xcorr(abs(refHRRProfile), abs(rangeProfile));
        % Find index of peak correlation value
        PeakIndex = find(CorrelationOut==max(CorrelationOut));
        % Calculate shift between ref peak and this range profile peak
        shift = PeakIndex - PeakIndexRef;
        shifts(indx)=shift;

        % Shift range profile to align with ref
        shiftedRangeProfile = circshift(rangeProfile,shift);
        % Add shifted profile to array
        RA_HRR_profiles(indx,:) = shiftedRangeProfile;
    end
end