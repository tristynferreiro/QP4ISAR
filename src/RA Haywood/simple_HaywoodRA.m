%% Implementation of the Haywood RA algorithm
% Initial explanation for understanding:
RangeProfiles = [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 
    0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0];
figure; imagesc(abs(RangeProfiles))
colorbar

%% Real implementation
% This uses the Haywood Range-Alignment algorithm (no autofocus). Haywood's
% algorithm begins with the simple correlation
[ShiftedRangeProfiles,shifts] = correlationRA(RangeProfiles);

%disp(ShiftedRangeProfiles)
figure; imagesc(abs(ShiftedRangeProfiles))
colorbar
