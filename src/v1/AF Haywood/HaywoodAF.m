function [AF_RA_HRRP] = HaywoodAF(RA_HRRP)
% Implements the Haywood autofocus algorithm
    % This approach is based on the dominant scatterer algorithm (DSA). 

[numProfiles,numRangeBins] = size(RA_HRRP);

% Step 1: Calculate Eq A.8 of Zyweck's appendix - amplitude variance
amplitudeVariance = (1/(numProfiles-1))*var(abs(RA_HRRP), [],1); % 1xM matrix

% Step 2: Identify scatterer (range bin) that satifies the dominant scatterer criteris. This is the minimum amplitude varying range bin in each range bin row
% Criteria 1: power of scatterer>average power - Eq A.10 of Zyweck's appendix
powerScatterer = sum(abs(RA_HRRP).^2,1);
averagePowerScatterer = mean(powerScatterer);
% find possible scatterers where power of scatterer > average power
candidateScatterersIdx = find(powerScatterer>averagePowerScatterer); % array of matches

% Criteria 2: candidate with minimum variance is dominant scatterer (DS) - Eq A.9 of Zyweck's appendix
[~,varianceDSidx] = min(amplitudeVariance(candidateScatterersIdx)); % returns index of match in candidateScatterersIdx
DSidx = candidateScatterersIdx(varianceDSidx); % range bin number of DS

% Step 3: Calculate phase differences between reference and other profiles - Eq A.11 of Zyweck's appendix
phaseHistoryDS = angle(RA_HRRP(:,DSidx)); % N x 1 matrix
disp(DSidx)
% phase_ref = abs(RA_HRRP(:,DSidx)); 
% phase_diff = abs(RA_HRRP) - phase_ref; % Gives Nx1 array

% Step 4: apply phase differences to each HRRP - Eq A.12 of Zyweck's appendix
complexConjugateDS = exp(-1i*phaseHistoryDS);
compensationMatrix = repmat(complexConjugateDS,1,numRangeBins);
AF_RA_HRRP = RA_HRRP.*compensationMatrix; % Eq A.13 of Zyweck's appendix
%AF_RA_HRRP = RA_HRRP.*complexConjugateDS; 
end