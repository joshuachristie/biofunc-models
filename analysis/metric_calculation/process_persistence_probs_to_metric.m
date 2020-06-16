% script that calculates metric for biological function
% imports persistence probabilities (calculated by the
% C++ models in /models/wright_fisher/ that are called by the python
% scripts in /models/run/)
% for all models but haploid_single_environment, it calls
% /analysis/PID/redAsMinIComponentTotal.m, which performs partial
% information decomposition.
% Processes entire set of persistence probabilities from
% /data/persistence_probabilities and prints data to
% /data/calculated_metric_values
addpath('../PID');
subdirs = dir('../../data/persistence_probabilities');
subdirs = subdirs(~ismember({subdirs.name}, {'.', '..'}));

for i = 1:numel(subdirs)
    if ~strcmp(subdirs(i).name, 'HSE') % HSE doesn't require PID
        path_to_folder = strcat(subdirs(i).folder, '/', subdirs(i).name);
        list_of_files = dir(path_to_folder);
        list_of_files = list_of_files(~ismember({list_of_files.name}, {'.', '..'}));
        
        % csv has four floats and contains persistence probabilities according to the following scheme (0-based indexing):
        
        % DSE
        % 0th element: neutral drift
        % 1st element: heterozygote effect alone
        % 2nd element: homozygote effect alone
        % 3rd element: focal situation (effect of both heterozygote and homozygote)
        
        % HTEOE
        % 0th element: neutral drift
        % 1st element: selection_coefficient_A2 alone
        % 2nd element: selection_coefficient_A1 alone
        % 3rd element: focal situation (effect of both selection_coefficient_A1 and selection_coefficient_A2)
        
        % HTE
        % 0th element: neutral drift
        % 1st element: selection_coefficient_A2 alone
        % 2nd element: selection_coefficient_A1 alone
        % 3rd element: focal situation (effect of both selection_coefficient_A_env_1 and selection_coefficient_A_env_2)
        
        % source 1 refers to the homozygote (DSE) and A1 (HTEOE/HTE)
        % source 2 refers to the heterozygote (DSE) and A2 (HTEOE/HTE)
        
        % calculate overall metric (non-decomposed)
        for j = 1:numel(list_of_files)
            persistence_probs = csvread(strcat(path_to_folder, '/', list_of_files(j).name));
            % calculate_metric takes prob due to drift and due to focal situation
            metric = calculate_metric(persistence_probs(1), persistence_probs(4));
            % get normalised weighted PID (R, U1, U2, C)
            pTable = [
                0, 0, 1, persistence_probs(1) / 4; 0, 0, 0, (1-persistence_probs(1)) / 4;...
                0, 1, 1, persistence_probs(2) / 4; 0, 1, 0, (1-persistence_probs(2)) / 4;...
                1, 0, 1, persistence_probs(3) / 4; 1, 0, 0, (1-persistence_probs(3)) / 4;...
                1, 1, 1, persistence_probs(4) / 4; 1, 1, 0, (1-persistence_probs(4)) / 4
                ];

            normalised_PID_both_sources_present = redAsMinIComponentTotal(pTable);
            % get decomposed metric and write to data/
            metric_decomposition = metric * normalised_PID_both_sources_present;
            csvwrite(sprintf('../../data/calculated_metric_values/%s/%s', subdirs(i).name, ...
                regexprep(list_of_files(j).name, '.csv', '_metric_decomposition.csv')), metric_decomposition);
        end
        
    else
        % HSE
        % csv has two floats (doesn't require PID)
        % 0th element: neutral drift
        % 1st element: selection_coefficient
        
        path_to_folder = strcat(subdirs(i).folder, '/', subdirs(i).name);
        list_of_files = dir(path_to_folder);
        list_of_files = list_of_files(~ismember({list_of_files.name}, {'.', '..'}));
        
        % calculate metric and write to data/
        for j = 1:numel(list_of_files)
            persistence_probs = csvread(strcat(path_to_folder, '/', list_of_files(j).name));
            metric = calculate_metric(persistence_probs(1), persistence_probs(2));
            csvwrite(sprintf('../../data/calculated_metric_values/%s/%s', subdirs(i).name, ...
                regexprep(list_of_files(j).name, '.csv', '_metric.csv')), metric);
        end
        
    end
end

