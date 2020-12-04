% Seth Talyansky, 2018-2020
% sethtal@stanford.edu
%*****************************************************
% For work stemming from use of this code, please cite
% Talyansky & Brinkman (2020) "Dysregulation of excitatory neural firing replicates physiological and functional changes in aging visual cortex", PLoS Computational Biology (2020).
% Code is adapted from RunVisionNetworkSimulation.m of E-I Net by Paul King
% (https://github.com/paulking/ei_net)
% E-I Net reference paper is King PD, Zylberberg J, DeWeese MR (2013). Inhibitory interneurons decorrelate
%        excitatory cells to drive sparse code formation in a spiking model of V1.
%        J Neuroscience.
%*****************************************************
%Inputs: model, grating stimulus patch-series matrix, number of
%orientations of grating stimuli, number of frames per orientation (angle)
%Purpose: Compute selectivity of the neurons in the network and store as an
%attribute of the model (m.selectivity)

function m = testSelectivity(m, ginputDataMaster, numOrientations, numfpa)
i_cg_input        = FindByName(m.cellGroup, 'input');
i_cg_output       = m.outputCellGroupId;
m.learningRate = 0;
m.numIterationsPerSample = 200; %longer for grating stimulus presentations; 50 during training
response = zeros(m.cg_V1e.numCells + m.cg_V1i.numCells, numOrientations); %spike counts in response to stimuli of different orientations
for theta = 0:numOrientations-1
    for f = 1:numfpa
        inputData = squeeze(ginputDataMaster(theta+1,f,:,:));
        m.snapshot.inputData = inputData;
        if f > 1
            initialState = finalState;
        else
            initialState = struct();
        end
        
        initialState.y{i_cg_input} = inputData * (m.inputScale * m.simTimeStep);
        
        % simulate the network in batch
        [m, finalState] = NetModel_UpdateFast(m, m.numIterationsPerSample, initialState);
        m.finalState = finalState;
        numSpikesE = mean(sum(finalState.y_history{i_cg_output}, 3), 2);
        numSpikesI = mean(sum(finalState.y_history{i_cg_output+1}, 3), 2);
        for s = 1:size(numSpikesE, 1)
            response(s, theta+1) = response(s, theta+1) + numSpikesE(s);
        end
        
        for s = 1:size(numSpikesI, 1)
            response(s+m.cg_V1e.numCells, theta+1) = response(s+m.cg_V1e.numCells, theta+1) + numSpikesI(s);
        end
       
    end
end
m.response = response;
m = orientationSelAnalysis(m); %compute selectivity based on response matrix
fprintf('mean selectivity of excitatory neurons is %d\n', mean(m.selectivity.exc)); 
end