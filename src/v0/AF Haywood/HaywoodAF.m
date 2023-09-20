function [AF_RA_HRRP] = HaywoodAF(RA_HRRP)
% Implements the Haywood autofocus algorithm
    % This approach is based on the dominant scatterer algorithm (DSA). 

    numRangeBins = size(RA_HRRP,2);
    
    %% Step 1: Calculate Eq A.8 of Zyweck's appendix - amplitude variance
    amplitudeVariance = var(abs(RA_HRRP), 1); % 1xM matrix
    
    %% Step 2: Identify scatterer (range bin) that satifies the dominant scatterer criteria.
    
    % Criteria 1: power of scatterer>average power - Eq A.10 of Zyweck's appendix
    powerScatterer = sum(abs(RA_HRRP).^2,1);
    averagePowerScatterer = mean(powerScatterer);
    % find possible scatterers where power of scatterer > average power
    scalingFactor = 1;
    candidateScatterersIdx = find(powerScatterer>scalingFactor*averagePowerScatterer); % array of matches
    
    % Criteria 2: candidate with minimum variance is dominant scatterer (DS) - Eq A.9 of Zyweck's appendix
    [~,varianceDSidx] = min(amplitudeVariance(candidateScatterersIdx)); % returns index of match in candidateScatterersIdx
    DSidx = candidateScatterersIdx(varianceDSidx); % range bin number of DS
    
    %% Step 3: Calculate phase differences - Eq A.11 of Zyweck's appendix
    phaseHistoryDS = angle(RA_HRRP(:,DSidx)); % N x 1 matrix
    
    %% Step 4: apply phase differences to each HRRP - Eq A.12 of Zyweck's appendix
    complexConjugateDS = exp(-1i*phaseHistoryDS);
    compensationMatrix = repmat(complexConjugateDS,1,numRangeBins);
    AF_RA_HRRP = RA_HRRP.*compensationMatrix; % Eq A.13 of Zyweck's appendix
end