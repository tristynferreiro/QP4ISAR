function [AF_RA_HRRP] = YuanAF_noise(RA_HRRP)
% Implements the Yuan autofocus algorithm
    % This approach is based on the dominant scatterer algorithm (DSA). 
    % This version considers all possible scatterers from all range bins
    % (even ones in which the target does not feature in).

    num_range_bins = size(RA_HRRP,2);
    
    %% Step 1: Calculate mean and variance
    amplitude_mean = mean(abs(RA_HRRP),1);
    amplitude_variance = var(abs(RA_HRRP),1);
    
    %% Step 2:  find candidate scatterers: 
    % Yuan's candidate scatterers with var/(var+mean) < 0.16
    criteria = amplitude_variance./(amplitude_variance+amplitude_mean.^2);
    candidate_scatterers_idx  = find(criteria< 0.16); % profile numbers
    
    %% Step 3:  Choose smallest 11 (preferred) but can choose number between 6-18
    numScatterers = 11; % ideally 11 otherwise value in range 6-18
    [~,candidate_scatterers_idx_min] = mink(criteria(candidate_scatterers_idx), numScatterers);
    DSidx = candidate_scatterers_idx(candidate_scatterers_idx_min); % get range bin numbers
    
    %% Step 4:  Determine constant phase shift for N pulses
    ref_bins = RA_HRRP(1,DSidx); % reference profile
    product_vector = conj(ref_bins).* RA_HRRP(:,DSidx);
    % we want the average phase difference across all range bins in each profile
    phase_shifts = angle(mean(product_vector,2)); 
    
    %% Step 5:  Apply phase shift
    compensation_angle = exp(-1i*phase_shifts);
    compensation_matrix = repmat(compensation_angle,1,num_range_bins);
    AF_RA_HRRP = RA_HRRP.*compensation_matrix;

end