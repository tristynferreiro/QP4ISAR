function [AF_RA_HRRP] = YuanAF_v0(RA_HRRP)
% Implements the Yuan autofocus algorithm
    % It is a type of dominant scatterer algorithm (DSA) that uses multiple 
    % dominant scatterres (DS) to achieve phase correction. Candidate 
    % scatterers are selected based on a threshold value of 0.16 
    % (step 2). From these candidate scatterers a number of DSs are
    % selected (step 3). The phases of the DSs are used to calculate a 
    % compensation factor use to correct all other scatterers (step 4).    
    %
    % Chao Yuan and David P. Casasent, https://doi.org/10.1117/1.1425789
    %
    % Version: v0

    %% Variabless
    num_range_bins = size(RA_HRRP,2);
    
    %% Step 1: Calculate mean and variance
    amplitude_mean = mean(abs(RA_HRRP),1);
    amplitude_variance = var(abs(RA_HRRP),1);
    
    %% Step 2:  find candidate scatterers:
    % Identify candidate scatterers using Yuan's threshold
    criteria = amplitude_variance./(amplitude_variance+amplitude_mean.^2);
    candidate_scatterers_idx  = find(criteria< 0.16); % profile numbers
    
    %% Step 3:  Choose smallest 11 (preferred) but can choose number between 6-18
    num_scatterers = 11; % ideally 11 otherwise value in range 6-18

    [~,candidate_scatterers_idx_min] = mink(criteria(candidate_scatterers_idx), num_scatterers);
    DS_idx = candidate_scatterers_idx(candidate_scatterers_idx_min); % get range bin numbers
    
    %% Step 4:  Determine constant phase shift for N pulses
    ref_bins = RA_HRRP(1,DS_idx); % reference profile
    product_vector = conj(ref_bins).* RA_HRRP(:,DS_idx);
    % calculate average phase difference across all range bins in each profile
    phase_shifts = angle(mean(product_vector,2)); 

    %% Step 5:  Apply phase shift
    compensation_angle = exp(-1i*phase_shifts);
    compensation_matrix = repmat(compensation_angle,1,num_range_bins);
    AF_RA_HRRP = RA_HRRP.*compensation_matrix;

end