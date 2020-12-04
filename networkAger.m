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
% Comments with "(King)" are preserved from original E-I Net code
%****************************************************


%Inputs: natural image matrix (image height x image width x number of
%images), grating image matrix (same format), parent directory within which to
%save network models at different ages
%Purpose: create the network model and simulate it through different ages
function networkAger(IMAGES, gratings, basedir)
format long e;
addpath('util');
addpath('neurosim');
numImages = size(IMAGES, 3);
gnumImages = size(gratings, 3);

%%
%this section is adapted from demo3_Train_EINetFast.m and
%RunVisionNetworkSimulation.m of E-I Net
%setting basic parameters of network
simParams = struct();
simParams.numInputSamples              =  400000; %5000 * number of training loops (80)
simParams.model.modelType              = 'V1eV1i';
simParams.model.learningRate = 0.4;
simParams.model.autoFigures            = {'STA_E', 'STA_I'};
simParams.model.stats.numSTASamples    = 5000;            % average 5000 samples when computing STA
simParams.model.inputDims = [8,8];
defaultValues.numInputSamples       = 10000;
defaultValues.printSummary          = false;       % print summary report at end
defaultValues.model.modelType       = [];         % no default model type
defaultValues.model.inputPreprocess.csFilter             = false;
defaultValues.model.autoAnneal.wait                      = 0;
defaultValues.model.autoAnneal.minLearningRate           = .01;
defaultValues.model.autoAnneal.tolerance                 = .001;
defaultValues.model.stats.measure.resErr_out.measureType = 'resError';
defaultValues.model.stats.measure.resErr_out.sourceName  = 'cg_output.in_input';
simParams = ApplyDefaultValues(simParams, defaultValues);
%%

numVariants = 5;
freqtests = 5; %the interval of loops at which the network selectivity is examined
youthlength = 30; %the length in loops of the maturation phase, before aging is initiated
etsrincrement = 0.002; %excitatory target spike rate increment per training loop (unit of age)
 
%variant 1 is normal (control) aging case, variant 2 is case with input
%weights frozen during aging (loops 31 through 80); variant 3 is case with
%lateral weights frozend during aging; variant 4 is case with both input
%and lateral weights frozen during aging; variant 5 is case with learning
%rate reduced 10x during aging and grating presentation time extended 10x
%(critical learning period theory test); variant 6 is case with learning
%rate reduced 10x during aging and excitatory target spike rate reduced
%10x(critical learning period theory test)
for i = 1:numVariants
    mkdir(sprintf('%s/model%d', basedir, i));    
end
%%
%this section is adapted from demo3_TrainEINetFast.m of E-I Net
models = cell(1, numVariants);
for i = 1:numVariants
    model = NetModel_InitV1eV1i(simParams.model);
    model.simNoReset = true; %set to false in default E-I Net
    model.learningRate = 0.4; %learning rate used in demo3_TrainEINetFast.m of E-I Net package
    model.stats.printStatsEvery = 0; %eliminates periodic reporting about training to console
    model.cellGroup{3}.inputBlock(1).connectionType = 'continuous'; %feed image input directly to inhibitory neurons, too; E-I Net default for this connection is 'disabled'
    model.cellGroup{2}.targetSpikeRate = 0.01; %E-I Net default is 0.05
    models{1, i} = model;
end
%%
%this section is adapted from RunVisionNetworkSimulation.m of E-I Net


zLength = 50; %number of frames presented in series per iteration of training
patchDims = model.inputDims;
numSamplesPerBatch = model.numSamplesPerBatch; %since all models have same number of samples per batch
numInputSamples = simParams.numInputSamples;
ind_r_c = zeros(numSamplesPerBatch, 3); %image indices, row indices, column indices
i_cg_input        = FindByName(model.cellGroup, 'input');
i_cg_output       = model.outputCellGroupId;

numLoops = numInputSamples/(zLength*numSamplesPerBatch);

numOrientations = 4;
numfpa = floor(gnumImages/numOrientations); %num frames per angle orientation
gind_r_c = zeros(numSamplesPerBatch, 2); %image indices, row indices, column indices


for l = 1:numSamplesPerBatch
    gind_r_c(l, 1)  = ceil( (size(gratings, 1)-patchDims(1)) * rand() );
    gind_r_c(l, 2)  = ceil( (size(gratings, 2)-patchDims(2)) * rand() );
end


ginputDataMaster = zeros(numOrientations, numfpa, prod(patchDims), numSamplesPerBatch);
for theta = 0:numOrientations-1
    inputData = zeros(prod(patchDims), numSamplesPerBatch);
    for f = 1:numfpa
        for j = 1:numSamplesPerBatch
            r = gind_r_c(j, 1);
            c = gind_r_c(j, 2);
            imagePatch = gratings(r:r+patchDims(1)-1, c:c+patchDims(2)-1, theta*numfpa+f);
            inputData(:,j) = imagePatch(:);
        end
        
        % normalize to zero mean and unit variance (King)
        inputData = bsxfun(@minus, inputData, mean(inputData));
        inputData = bsxfun(@times, inputData, 1./std(inputData));
        ginputDataMaster(theta+1,f,:,:) = inputData;
    end
end

for n = 1:numLoops   
    if n > youthlength
        for v = 1:numVariants
            model = models{1, v};            
                if n == youthlength + 1
                    if v == 2
                        model.cellGroup{2}.inputBlock(1).learningRate = 0;
                        model.cellGroup{3}.inputBlock(1).learningRate = 0;
                    elseif v == 3
                        model.cellGroup{2}.inputBlock(3).learningRate = 0;
                        model.cellGroup{3}.inputBlock(3).learningRate = 0;
                        model.cellGroup{3}.inputBlock(2).learningRate = 0;
                    elseif v == 4
                        model.cellGroup{2}.inputBlock(1).learningRate = 0;
                        model.cellGroup{3}.inputBlock(1).learningRate = 0;
                        model.cellGroup{2}.inputBlock(3).learningRate = 0;
                        model.cellGroup{3}.inputBlock(3).learningRate = 0;
                        model.cellGroup{3}.inputBlock(2).learningRate = 0;
                    elseif v == 5
                        model.learningRate = model.learningRate/10;
                        model.numIterationsPerSample = model.numIterationsPerSample*10;
                    end
                end
            model.cellGroup{2}.targetSpikeRate = model.cellGroup{2}.targetSpikeRate + etsrincrement;
            models{1, v} = model;
        end
    end
    fprintf('Beginning simulation step of loop %d of %d\n', n, numLoops);
    for l = 1:numSamplesPerBatch %selection of sample (patch) coordinates
        ind_r_c(l, 1)  = ceil( (floor(numImages/2) - zLength)      * rand() );
        ind_r_c(l, 2)  = ceil( (size(IMAGES, 1)-patchDims(1)) * rand() );
        ind_r_c(l, 3)  = ceil( (size(IMAGES, 2)-patchDims(2)) * rand() );
    end
    for s = 1:zLength
        inputData = zeros(prod(patchDims), numSamplesPerBatch);
        for j = 1:numSamplesPerBatch
            r = ind_r_c(j, 2);
            c = ind_r_c(j, 3);
            imagePatch = IMAGES(r:r+patchDims(1)-1, c:c+patchDims(2)-1, ind_r_c(j, 1)+s-1);
            inputData(:,j) = imagePatch(:);
        end
        % normalize to zero mean and unit variance (King)
        inputData = bsxfun(@minus, inputData, mean(inputData));
        inputData = bsxfun(@times, inputData, 1 ./ std(inputData));
        for v = 1:numVariants
            model = models{1, v};
            model.snapshot.inputData = inputData;
            
            if isfield(model, 'finalState')
                % special case simNoReset mode: copy state if after the
                % first iteration (King)
                initialState = model.finalState;
            else
                initialState = struct();
            end
            
            initialState.y{i_cg_input} = inputData * (model.inputScale * model.simTimeStep);
            
            % simulate the network in batch (King)
            [model, finalState] = NetModel_UpdateFast(model, model.numIterationsPerSample, initialState); %NetModel_UpdateFast is training with current; NetModel_UpdateFastControl is training without current
            model.finalState = finalState;
            models{1, v} = model;
        end
    end
    if mod(n, freqtests) == 0
        fprintf('Beginning testing/data collection step of loop %d of %d\n', n, numLoops);
        for v = 1:numVariants
            %orientation selectivity testing
            model = models{1, v};
            m = model; %duplicating model such that selectivity testing does not alter original object
            m.learningRate = 0;            
            m.numIterationsPerSample = 200; %longer for grating stimulus presentations           
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
                    
                    % simulate the network in batch (King)
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
            model.response = m.response;
            m = orientationSelAnalysis(m); %compute selectivity based on response matrix
            model.selectivity = m.selectivity;
            save(sprintf('%s/model%d/ln_%d.mat', basedir, v, n), 'model');
        end
        
    end
    
end
%%
end
