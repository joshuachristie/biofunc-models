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

subdirs = dir('../../data/persistence_probabilities');
subdirs = subdirs(~ismember({subdirs.name}, {'.', '..'}));

for i = 1:numel(subdirs)
    if ~strcmp(subdirs(i).name, 'HSE') % HSE doesn't require PID
        path_to_folder = strcat(subdirs(i).folder, '/', subdirs(i).name);
        list_of_files = dir(path_to_folder);
        list_of_files = list_of_files(~ismember({list_of_files.name}, {'.', '..'}));
        
        for j = 1:numel(list_of_files)            
            persistence_probs = csvread(strcat(path_to_folder, '/', list_of_files(j).name)); 
        end
            
        % csv has four floats and contains persistence probabilities according to the following scheme:
        
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

        pTable = [
            0, 0, 1, persistence_probs(1) / 4; 0, 0, 0, (1-persistence_probs(1)) / 4;...
            0, 1, 1, persistence_probs(2) / 4; 0, 1, 0, (1-persistence_probs(2)) / 4;...
            1, 0, 1, persistence_probs(3) / 4; 1, 0, 0, (1-persistence_probs(3)) / 4;...
            1, 1, 1, persistence_probs(4) / 4; 1, 1, 0, (1-persistence_probs(4)) / 4
            ];
        
        % call redAsMinIComponentTotal(pTable) to get normalised PID scores
       % call function to calculate metric
       % weight metric by PID
    else
        % call function to calculate metric for HSE (no PID)
    end
end
        