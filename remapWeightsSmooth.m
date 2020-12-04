% Seth Talyansky, 2018-2020
% sethtal@stanford.edu
%*****************************************************
% For work stemming from use of this code, please cite
% Talyansky & Brinkman, "Dysregulation of excitatory neural firing replicates physiological and functional changes in aging visual cortex", PLoS Computational Biology (2020).
%*****************************************************
%Inputs: young network model and aged network model
%Purpose: produce network model with input weights that have old magnitudes
%but young spatial distribution (implements Eq. (12) of Talyansky &
%Brinkman (2020))
%To produce network model with input weights that have young magnitudes but
%old spatial distribution simply reorder inputs
%Curve Fitting Toolbox must be installed

function young_aged = remapWeightsSmooth(young, old)
young_aged = old;
for cg = 2:3 %excitatory input weights and then inhibitory input weights
    young_input_weights = young.cellGroup{cg}.inputBlock(1).weight(:);
    young_cdf_discrete = cdfplot(young_input_weights);
    young_cdf_x = young_cdf_discrete.XData';
    infinds = find(abs(young_cdf_x) == Inf);
    young_cdf_x(infinds) = [];
    young_cdf_y = young_cdf_discrete.YData';
    young_cdf_y(infinds) = [];
    young_cdf = fit(young_cdf_x, young_cdf_y, 'smoothingspline'); %cumulative distribution function of young input weights
    
    old_input_weights = old.cellGroup{cg}.inputBlock(1).weight(:);
    old_cdf_discrete = cdfplot(old_input_weights);
    old_cdf_x = old_cdf_discrete.XData';
    infinds = find(abs(old_cdf_x) == Inf);
    old_cdf_y = old_cdf_discrete.YData';
    old_cdf_x(infinds) = [];
    old_cdf_y(infinds) = [];
    old_cdf_inv = fit(old_cdf_y, old_cdf_x, 'smoothingspline'); %inverse cumulative distribution function of old input weights
    
    for i = 1:numel(young_input_weights)
        young_aged.cellGroup{cg}.inputBlock(1).weight(i) = old_cdf_inv(young_cdf(young_input_weights(i)));
    end
end
end