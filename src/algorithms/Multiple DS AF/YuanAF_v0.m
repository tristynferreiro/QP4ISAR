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
    %
    % Revision:
    %   - Intemediary step to select the reference profile. This was not
    %       suggested in the paper - it is a new addition.

    %% Variabless
    num_range_bins = size(RA_HRRP,2);
    
    %% Step 1: Calculate mean and variance
    amplitude_mean = mean(abs(RA_HRRP),1);
    amplitude_variance = var(abs(RA_HRRP),1);
    
    %% Step 2:  find candidate scatterers:
    % Identify candidate scatterers using Yuan's threshold
    dispersion = amplitude_variance./(amplitude_variance+amplitude_mean.^2);
    candidate_scatterers_idx  = find(dispersion< 0.16); % profile numbers

    %% Step 3:  Choose smallest 11 (preferred) but can choose number between 6-18
    num_scatterers = 11; % ideally 11 otherwise value in range 6-18

    [~,candidate_scatterers_idx_min] = mink(dispersion(candidate_scatterers_idx), num_scatterers);
    DS_idx = candidate_scatterers_idx(candidate_scatterers_idx_min); % get range bin numbers

    % Plot the DS selection to validate the selection is correct
    figure; plot(dispersion,'-*'); hold on;
    yline(0.16,'-g');
    plot(candidate_scatterers_idx, dispersion(candidate_scatterers_idx), 'ok', 'MarkerSize', 8);
    plot(DS_idx, dispersion(DS_idx), 'or', 'MarkerSize', 15);
    xlabel('Range Bin Number'); ylabel('Scatterer Dispersion');
    title('Candidate DS selection.');
    legend('Scatterer Dispersion', 'Threshold','Candidate scatterers','Dominant Scatterers');
    hold off;
    matlab2tikz() % Only uncomment to Save the figure as LaTeX compatible plot
    
    %% Intemediary step
    % This was an addition to the original algorithm described by Yuan.
    % Choose the reference profile as the profile where the sum of all DS
    % power is the greatest.
    [~,ref_HRRP_num] = max(sum(abs(RA_HRRP(:,DS_idx)).^2,2));
    
    % % Plot the selected reference profile
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

    % Angle before phase correction
    % scatterer =11;
    % figure; plot(angle(RA_HRRP(:,DS_idx(scatterer))));
    % xlabel('Profile Number'); ylabel('Phase Angle');
    % title(['Dominant Scatterer ',num2str(scatterer), ' Phase']);
    %matlab2tikz() % Only uncomment to Save the figure as LaTeX compatible plot
    
    % Angle after phase correction
    % figure; plot(angle(AF_RA_HRRP(:,DS_idx(scatterer))));
    % xlabel('Profile Number'); ylabel('Phase Angle');
    % title(['Dominant Scatterer ',num2str(scatterer), ' Corrected Phase']);
    %matlab2tikz() % Only uncomment to Save the figure as LaTeX compatible plot

end