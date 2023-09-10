%% Implementation of the Haywood RA algorithm
RangeProfiles = [1+0.3i 0 0 0 0 0 0 0; 0 1-0.1i 0 0 0 0 0 0; 
    0 0 1+0.7i 0 0 0 0 0; 0 0 0 1-1.2i 0 0 0 0; 0 0 0 0 1-1.2i 0 0 0;
    0 0 0 0 0 1-1.2i 0 0; 0 0 0 0 0 0 1+1.2i 0; 0 0 0 0 0 0 0 1-1.2i];
figure; imagesc(abs(RangeProfiles))
colorbar
%% Real implementation of RA
% Range Align the HRR profiles using Haywood method
ref_profile_number = 1;
RA_HRRP = HaywoodRA(RangeProfiles, ref_profile_number);

figure; imagesc(abs(RA_HRRP))
colorbar
%% Real Haywood AF implementation
[numProfiles,numRangeBins] = size(RangeProfiles);
% Step 1: Calculate Eq A.8 of Zyweck's appendix - amplitude variance
amplitudeVariance = 1/(numProfiles-1)*var(abs(RA_HRRP), [],1); % 1xM matrix

% Step 2: Identify scatterer (range bin) that satifies the dominant scatterer criteris. This is the minimum amplitude varying range bin in each range bin row
% Criteria 1: power of scatterer>average power - Eq A.10 of Zyweck's appendix
powerScatterer = sum(abs(RA_HRRP).^2,1);
averagePowerScatterer = mean(powerScatterer);
% find possible scatterers where power of scatterer > average power
candidateScatterersIdx = find(powerScatterer>averagePowerScatterer); % array of matches

% Criteria 2: minimum variance is dominant scatterer (DS) - Eq A.9 of Zyweck's appendix
[varianceDS,varianceDSidx] = min(amplitudeVariance(candidateScatterersIdx)); % returns index

% Step 3: Calculate phase differences between reference and other profiles - Eq A.11 of Zyweck's appendix
phaseHistoryDS = angle(RA_HRRP(:,varianceDSidx)); % N x 1 matrix

% Step 4: apply phase differences to each HRRP - Eq A.12 of Zyweck's appendix
complexConjugateDS = exp(-1i*phaseHistoryDS);
% compensationMatrix = repmat(complexConjugateDS,1,numRangeBins)
% AF_RA_HRRP_1 = RA_HRRP.*compensationMatrix
AF_RA_HRRP = RA_HRRP.*complexConjugateDS; % Eq A.13 of Zyweck's appendix

figure; imagesc(abs(AF_RA_HRRP))
colorbar


