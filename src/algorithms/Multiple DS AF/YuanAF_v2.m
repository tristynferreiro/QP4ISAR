function [AF_RA_HRRP] = YuanAF_v2(RA_HRRP)
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
    % Version: v2
    %
    % Revisions: 
    %   - Intemediary step to select the reference profile. This was not
    %       suggested in the paper - it is a new addition.
    %   - Added if statement to set number of DSs to highest number
    %     possible in cases where a smaller number of candidate scatterers
    %     are available.
    %   - NEW: A threshold was introduced to ensure that only scatterers 
    %       within the object are considered. i.e. no noisy scatterers are 
    %       used. This prevents false triggering in DS selection 
    
    %% Variables
    num_range_bins = size(RA_HRRP,2);
    
    %% Step 1: Calculate mean and variance
    amplitude_mean = mean(abs(RA_HRRP),1);
    amplitude_variance = var(abs(RA_HRRP),1);
    
    %% Intemediary step: Scatterer filtering
    % Threshold to get object range bins
    no_noise_scatterers_idx = find(amplitude_mean>0.5.*mean(amplitude_mean));
    no_noise_scatterers_idx = min(no_noise_scatterers_idx):max(no_noise_scatterers_idx);

    % Plot noise filtered
    % figure; plot(amplitude_mean,'-*'); hold on;
    % yline(0.5.*mean(amplitude_mean),'-g');
    % plot(no_noise_scatterers_idx,amplitude_mean(no_noise_scatterers_idx),'or');
    % xlabel('Range Bin Number'); ylabel('Scatterer Amplitudes');
    % title('Filterered scatterers');
    % legend('Scatterer Amplitudes', 'Mean of Scatterer Amplitudes','Object Scatterers');
    % hold off;
    % matlab2tikz() % Only uncomment to Save the figure as LaTeX compatible plot

    %% Step 2:  find candidate scatterers: 
    % Identify candidate scatterers using Yuan's threshold
    dispersion = amplitude_variance./(amplitude_variance+amplitude_mean.^2);
    dispersion_thresh = dispersion(no_noise_scatterers_idx)< 0.16;
    candidate_scatterers_idx = no_noise_scatterers_idx(dispersion_thresh);

    %% Step 3:  Choose smallest x scatterers
    if(size(candidate_scatterers_idx,2)>11)
        num_scatterers = 11; % ideally 11 otherwise value in range 6-18
    else
        % set the number of scatterers as largest value value possible
        % given the number of candidate scatterers
        num_scatterers = size(candidate_scatterers_idx,2);
    end


    [~,candidate_scatterers_idx_min] = mink(dispersion(candidate_scatterers_idx), num_scatterers);
    DS_idx = candidate_scatterers_idx(candidate_scatterers_idx_min); % get range bin numbers

    % % Plot the DS selection to validate the selection is correct
    % figure; plot(dispersion,'-*'); hold on;
    % yline(0.16,'-g');
    % plot(candidate_scatterers_idx, dispersion(candidate_scatterers_idx), 'ok', 'MarkerSize', 8);
    % plot(DS_idx, dispersion(DS_idx), 'or', 'MarkerSize', 15);
    % xline([min(no_noise_scatterers_idx), max(no_noise_scatterers_idx)], '-m');
    % xlabel('Range Bin Number'); ylabel('Scatterer Dispersion');
    % title('Candidate DS selection.');
    % legend('Scatterer Dispersion','Threshold','Candidate scatterers','Dominant Scatterers', 'Object Scatterer Boundary');
    % hold off;
    % matlab2tikz() % Only uncomment to Save the figure as LaTeX compatible plot
    
    %% Intemediary step
    % This was an addition to the original algorithm described by Yuan.
    % Choose the reference profile as the profile where the sum of all DS
    % power is the greatest.
    [~,ref_HRRP_num] = max(sum(abs(RA_HRRP(:,DS_idx)).^2,2));
    
    % Plot the selected reference profile
    % customPurple = [0.5, 0, 0.5];
    % figure; plot(sum(abs(RA_HRRP(:,DS_idx)).^2,2),'-*'); hold on; 
    % plot(ref_HRRP_num,max(sum(abs(RA_HRRP(:,DS_idx)).^2,2)),'s', 'MarkerSize', 15,'Color', customPurple);
    % xline(ref_HRRP_num);
    % xlabel('Profile Number'); ylabel('Profile Power');
    % title('Reference Profile Selection');
    % legend('Sum of DS Power in each profile','Selected Profile');
    % grid on;  % Add grid lines
    % grid minor;  % Add minor grid lines
    % hold off;
    % matlab2tikz() % Only uncomment to Save the figure as LaTeX compatible plot
    
    %% Step 4:  Determine constant phase shift for N pulses
    ref_bins = RA_HRRP(ref_HRRP_num,DS_idx); % reference profile
    product_vector = conj(ref_bins).* RA_HRRP(:,DS_idx);
    % calculate average phase difference across all range bins in each profile
    phase_shifts = angle(mean(product_vector,2)); 

    %% Step 5:  Apply phase shift
    compensation_angle = exp(-1i*phase_shifts);
    compensation_matrix = repmat(compensation_angle,1,num_range_bins);
    AF_RA_HRRP = RA_HRRP.*compensation_matrix;

end