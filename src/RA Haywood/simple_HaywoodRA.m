%% Implementation of the Correlation RA algorithm
% Initial explanation for understanding:
RangeProfiles = [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 
    0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0];
figure; imagesc(abs(RangeProfiles))
colorbar

%RangeProfile1 = RangeProfiles(1,:);
%RangeProfile2 = RangeProfiles(2,:);
% THINK: want find a way to align 2 to 1. Intuitively we know to shift back
% by -1:
%ShiftedRangeProfile2 = circshift(RangeProfile2,-1);

%% Real implementation
% To find the shift value, we use the correlation between a reference
% profile and the profile we want to shift. The xcorr function returns an
% array of correlation values. The peak value will be at a different array
% index depending on the correlation between a profile and the reference.
% The location of the peak, when compared to the location of the peak in
% the autocorrelated return, indicates the range shift of the profile from
% the reference profile. The change in the postition of the peak in the 
% correlation of different profiles is because of the movement of the 
% target's dominant scatteres across range bins. If we compare the position 
% of the peak in the correlation with the autocorrelation then we can find 
% the range bin shift.

[ShiftedRangeProfiles,shifts] = correlationRA(RangeProfiles);

%disp(ShiftedRangeProfiles)
figure; imagesc(abs(ShiftedRangeProfiles))
colorbar
