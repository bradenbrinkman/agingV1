% Seth Talyansky, 2018-2020
% sethtal@stanford.edu
%*****************************************************
% For work stemming from use of this code, please cite
% Talyansky & Brinkman (2020) "Dysregulation of excitatory neural firing replicates physiological and functional changes in aging visual cortex", PLoS Computational Biology (2020).
%*****************************************************
%Inputs: model
%Purpose: Compute neuron selectivity measures from spike counts in response
%to oriented grating stimuli

function model = orientationSelAnalysis(model) %orientation selectivity analysis
%Formula for selectivity from Hua et al., "Functional
%degradation of visual cortical cells in old cats," Neurobiology of Aging
%(2006); Levick and Thibos, "Analysis of orientation bias in cat retina,"
%The Journal of Physiology (1982); Mazurek et al., "Robust quantification
%of orientation selectivity and direciton selectivity," Frontiers in Neural
%Circuits (2014)
response = model.response;
model.selectivity = struct();

selectivities = zeros(size(response, 1), 1);
for i = 1:size(response, 1)
    if sum(response(i, :)) == 0
        selectivities(i) = 0;
    else
        netx = 0;
        nety = 0;
        for j = 1:size(response, 2)
            netx = netx + cos((j-1)*pi/2)*response(i, j);
            nety = nety + sin((j-1)*pi/2)*response(i, j);
        end
        selectivities(i) = sqrt(netx^2 + nety^2)/sum(response(i, :));
    end
end

model.selectivity.exc = sort(selectivities(1:model.cg_V1e.numCells)); %selectivity indices of excitatory neurons
model.selectivity.inh = sort(selectivities((model.cg_V1e.numCells+1):numel(selectivities))); %selectivity indices of inhibitory neurons
end