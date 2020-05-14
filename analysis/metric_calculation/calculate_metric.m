% script that calculates metric for biological function
% imports persistence probabilities (calculated by the
% C++ models in /models/wright_fisher/ that are called by the python
% scripts in /models/run/)
% for all models but haploid_single_environment, it calls
% /analysis/PID/redAsMinIComponentTotal.m, which performs partial
% information decomposition.
% Prints data to /data/calculated_metric_values