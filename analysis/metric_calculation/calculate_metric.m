% function that returns the function metric value
% for a given neutral persistence probability
% pp_drift is the persistence probability for neutral drift
% pp_focal_case is the persistence probability for the case for which we
% are calculating the metric

% The metric is a real number in the range [-1, 1].
% We (axiomatically) define the following three points.
% metric = [-1, 0, 1]; persistence_probabilities = [0, ppd, 1]
% where ppd is the persistence_probability of pure drift
% (i.e. selection coefficient = 0)
% the metric, as implemented here, is simply a linear interpolation between these three points
% (note that one could an arbitrary non-decreasing function to join
% the points but as a default I'm focusing on the simplest case of a linear function)

% pp = persistence_probability

function metric = calculate_metric(pp_drift, pp_focal_case)
axiomatic_metric = [-1, 0, 1];
axiomatic_pp = [0, pp_drift, 1];
metric = interp1(axiomatic_metric, axiomatic_pp, pp_focal_case);
end