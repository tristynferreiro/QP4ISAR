function [AF_RA_HRRP] = YuanAF(RA_HRRP)
% Implements the Yuan autofocus algorithm
    % This approach is based on the dominant scatterer algorithm (DSA). 

numRangeBins = size(RA_HRRP,2);

% Step 1: Calculate mean and variance
amplitudeMean = mean(abs(RA_HRRP),1);
amplitudeVariance = var(abs(RA_HRRP),1);

% Step 2:  find candidate scatterers: 
% Power threshold to remove noise
powerScatterer = sum(abs(RA_HRRP).^2,1);
averagePowerScatterer = mean(powerScatterer);
% find possible scatterers where power of scatterer > average power
scalingFactor = 1;
noNoiseScatterers = find(powerScatterer>scalingFactor*averagePowerScatterer); % array of matches

% Yuan's candidate scatterers with var/(var+mean) < 0.16
criteria = amplitudeVariance./(amplitudeVariance+amplitudeMean.^2);
noNoiseCriteria=criteria(noNoiseScatterers);
idx  = find(noNoiseCriteria< 0.16); % profile numbers
candidateScatterersIdx = noNoiseScatterers(idx);

% Step 3:  Choose smallest 11 (preferred) but can choose number between 6-18
numScatterers = 11; % ideally 11 otherwise value in range 6-18
[~,candidateScatterersIdx_min] = mink(criteria(candidateScatterersIdx), numScatterers);
DSidx = candidateScatterersIdx(candidateScatterersIdx_min); % get range bin numbers

% Step 4:  Determine constant phase shift for N pulses
refBins = RA_HRRP(1,DSidx); % reference profile
product_vector = conj(refBins).* RA_HRRP(:,DSidx);
% we want the average phase difference across all range bins in each profile
phaseShifts = angle(mean(product_vector,2)); 

% Step 5:  Apply phase shift
compensationAngle = exp(-1i*phaseShifts);
compensationMatrix = repmat(compensationAngle,1,numRangeBins);
AF_RA_HRRP = RA_HRRP.*compensationMatrix;

end