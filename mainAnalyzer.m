% Seth Talyansky, 2018-2020
% sethtal@stanford.edu
%*****************************************************
% For work stemming from use of this code, please cite
% Talyansky & Brinkman (2020) "Dysregulation of excitatory neural firing replicates physiological and functional changes in aging visual cortex", PLoS Computational Biology (2020).
% Parts of the code are adapted from RunVisionNetworkSimulation.m of E-I Net by Paul King
% (https://github.com/paulking/ei_net)
% E-I Net reference paper is King PD, Zylberberg J, DeWeese MR (2013). Inhibitory interneurons decorrelate
%        excitatory cells to drive sparse code formation in a spiking model of V1.
%        J Neuroscience.
% Comments with "(King)" are preserved from original E-I Net code
% Curve Fitting (for remapWeightsSmooth.m) and Optimization (for receptive field Gabor character analysis at the end) Toolboxes must be installed
% Note: the violin plots provided in Talyansky & Brinkman (2020) are generated with the Python scripts included in this package; distributionPlot is included only for the convenience of generating similar plots directly in MATLAB

%****************************************************

%Driver function that produces image input if necessary, creates age-series of
%models, and collects statistics on those models
%naturalmovie is file containing cat camera film downloaded from https://zenodo.org/record/46481 (if file is not yet created, then provide path to folder of freshly downloaded images and make sure to include '/' at the end); resultsdir is destination for all model matrices and figures; if grating stimuli are not yet created, then instead of file name enter 'create'
%Usage e.g. mainAnalyzer('cat_movie_1.mat', 'gratings.mat', 'cm1/'); mainAnalyzer('cat_movie_7.mat', 'gratings.mat', 'cm7/'); mainAnalyzer('cat_movie_16.mat', 'gratings.mat', 'cm16/');
%or mainAnalyzer('catmovie1/', 'create', 'cm1'/) if whitened natural image
%matrix and grating stimulus image matrix do not yet exist
function mainAnalyzer(naturalmovie, gratingstimuli, resultsdir) 
format long e;
addpath util;
addpath neurosim;
addpath export_fig_container;
addpath distributionPlot;
addpath VisualStimulusToolbox-master/VisualStimulus;
randreset();
if ~exist(resultsdir, 'dir')
    mkdir(resultsdir);
end
%%
%this section adapted from Olshausen BA, Field DJ (1997)  Sparse Coding
%with an Overcomplete Basis Set: A Strategy Employed by V1?  Vision
%Research, 37: 3311-3325. Download: http://www.rctn.org/bruno/sparsenet/sparsenet.tar.gz
%STEP 1: produce whitened movie matrices

frameMargin = 4;

if isfolder(naturalmovie) %if movie file is folder of images and not yet whitened matrix
    
    listdir = dir(naturalmovie);
    M=numel(listdir)-2; %first two elements in directory are always '.' and '..'
    fn = listdir(3).name;
    image = imread(strcat(naturalmovie, '/', fn));
    x = size(image, 2); %image width in pixels; to be number of rows in matrix
    y = size(image, 1); %image height in pixels; to be number of columns in matrix
    
    IMAGES = zeros(y, x, M);
    
    [fx, fy]=meshgrid(-x/2:x/2-1,-y/2:y/2-1);
    rho=sqrt(fx.*fx/x^2+fy.*fy/y^2); %extension of initial code tailored to NxN images to account for possibility of dimensions of different sizes
    f_0=0.4;
    filt=rho.*exp(-(rho/f_0).^4);
    
    for i = 1:M
        fprintf('Whitening image %d of %d\n', i, M)
        fn = listdir(i+2).name;
        image = imread(strcat(naturalmovie, fn));
        If=fft2(image);
        IMAGES(:, :, i)=real(ifft2(If.*fftshift(filt)));
    end
    
    IMAGES=IMAGES/sqrt(mean(var(reshape(IMAGES, x*y, M))));
    
    % crop away a 4-pixel margin at the frame edge which contains FFT-induced
    % error (King)
    IMAGES = IMAGES(1+frameMargin:end-frameMargin, 1+frameMargin:end-frameMargin, :);
    save('cat_movie_1.mat', 'IMAGES'); %change number depending on movie
else
    imgs = load(naturalmovie);
    IMAGES = imgs.IMAGES;
end


if strcmp(gratingstimuli, 'create')
    y = size(IMAGES, 1) + 2 * frameMargin;
    x = size(IMAGES, 2) + 2 * frameMargin;
    grate = GratingStim([y x]);
    for orientation = (0:3)*45
        grate.add(10, orientation, [0.1 0.1], 0.1);
    end
    
    M=grate.length;
    gratings = zeros(y, x, M);
    
    [fx, fy]=meshgrid(-x/2:x/2-1,-y/2:y/2-1);
    rho=sqrt(fx.*fx/x^2+fy.*fy/y^2); %extension of original code tailored to NxN images to accommodate images with different dimensions
    f_0=0.4;
    filt=rho.*exp(-(rho/f_0).^4);
    
    grate_matrix = squeeze(grate.stim);
    for i = 1:M
        image = grate_matrix(:, :, i);
        If=fft2(image);
        gratings(:, :, i)=real(ifft2(If.*fftshift(filt)));
    end
    
    gratings=gratings/sqrt(mean(var(reshape(gratings, x*y, M))));
    gratings = gratings(1+frameMargin:end-frameMargin, 1+frameMargin:end-frameMargin, :);
    save('gratings.mat', 'gratings');
    
else
    grates = load(gratingstimuli);
    gratings = grates.gratings;
end
%%

%STEP 2: create and age model

networkAger(IMAGES, gratings, resultsdir);


%STEP 3: examine model

colors = parula(12);
colors(1:12,1:3) = colors([11:-1:1, 12],1:3);
loop_number = 30:5:80;
num_variants = 5;

%load model results

a30 = load(strcat(resultsdir, 'model1/ln_30.mat'));
a30 = a30.model;
a35 = load(strcat(resultsdir, 'model1/ln_35.mat'));
a35 = a35.model;
a40 = load(strcat(resultsdir, 'model1/ln_40.mat'));
a40 = a40.model;
a45 = load(strcat(resultsdir, 'model1/ln_45.mat'));
a45 = a45.model;
a50 = load(strcat(resultsdir, 'model1/ln_50.mat'));
a50 = a50.model;
a55 = load(strcat(resultsdir, 'model1/ln_55.mat'));
a55 = a55.model;
a60 = load(strcat(resultsdir, 'model1/ln_60.mat'));
a60 = a60.model;
a65 = load(strcat(resultsdir, 'model1/ln_65.mat'));
a65 = a65.model;
a70 = load(strcat(resultsdir, 'model1/ln_70.mat'));
a70 = a70.model;
a75 = load(strcat(resultsdir, 'model1/ln_75.mat'));
a75 = a75.model;
a80 = load(strcat(resultsdir, 'model1/ln_80.mat'));
a80 = a80.model;
control = [mean(a30.selectivity.exc), mean(a35.selectivity.exc), mean(a40.selectivity.exc), mean(a45.selectivity.exc), mean(a50.selectivity.exc), mean(a55.selectivity.exc), mean(a60.selectivity.exc), mean(a65.selectivity.exc), mean(a70.selectivity.exc), mean(a75.selectivity.exc), mean(a80.selectivity.exc)]; %mean orientation selectivity index of excitatory neurons at different ages 

i30 = load(strcat(resultsdir, 'model2/ln_30.mat'));
i30 = i30.model;
i35 = load(strcat(resultsdir, 'model2/ln_35.mat'));
i35 = i35.model;
i40 = load(strcat(resultsdir, 'model2/ln_40.mat'));
i40 = i40.model;
i45 = load(strcat(resultsdir, 'model2/ln_45.mat'));
i45 = i45.model;
i50 = load(strcat(resultsdir, 'model2/ln_50.mat'));
i50 = i50.model;
i55 = load(strcat(resultsdir, 'model2/ln_55.mat'));
i55 = i55.model;
i60 = load(strcat(resultsdir, 'model2/ln_60.mat'));
i60 = i60.model;
i65 = load(strcat(resultsdir, 'model2/ln_65.mat'));
i65 = i65.model;
i70 = load(strcat(resultsdir, 'model2/ln_70.mat'));
i70 = i70.model;
i75 = load(strcat(resultsdir, 'model2/ln_75.mat'));
i75 = i75.model;
i80 = load(strcat(resultsdir, 'model2/ln_80.mat'));
i80 = i80.model;

l30 = load(strcat(resultsdir, 'model3/ln_30.mat'));
l30 = l30.model;
l35 = load(strcat(resultsdir, 'model3/ln_35.mat'));
l35 = l35.model;
l40 = load(strcat(resultsdir, 'model3/ln_40.mat'));
l40 = l40.model;
l45 = load(strcat(resultsdir, 'model3/ln_45.mat'));
l45 = l45.model;
l50 = load(strcat(resultsdir, 'model3/ln_50.mat'));
l50 = l50.model;
l55 = load(strcat(resultsdir, 'model3/ln_55.mat'));
l55 = l55.model;
l60 = load(strcat(resultsdir, 'model3/ln_60.mat'));
l60 = l60.model;
l65 = load(strcat(resultsdir, 'model3/ln_65.mat'));
l65 = l65.model;
l70 = load(strcat(resultsdir, 'model3/ln_70.mat'));
l70 = l70.model;
l75 = load(strcat(resultsdir, 'model3/ln_75.mat'));
l75 = l75.model;
l80 = load(strcat(resultsdir, 'model3/ln_80.mat'));
l80 = l80.model;

il30 = load(strcat(resultsdir, 'model4/ln_30.mat'));
il30 = il30.model;
il35 = load(strcat(resultsdir, 'model4/ln_35.mat'));
il35 = il35.model;
il40 = load(strcat(resultsdir, 'model4/ln_40.mat'));
il40 = il40.model;
il45 = load(strcat(resultsdir, 'model4/ln_45.mat'));
il45 = il45.model;
il50 = load(strcat(resultsdir, 'model4/ln_50.mat'));
il50 = il50.model;
il55 = load(strcat(resultsdir, 'model4/ln_55.mat'));
il55 = il55.model;
il60 = load(strcat(resultsdir, 'model4/ln_60.mat'));
il60 = il60.model;
il65 = load(strcat(resultsdir, 'model4/ln_65.mat'));
il65 = il65.model;
il70 = load(strcat(resultsdir, 'model4/ln_70.mat'));
il70 = il70.model;
il75 = load(strcat(resultsdir, 'model4/ln_75.mat'));
il75 = il75.model;
il80 = load(strcat(resultsdir, 'model4/ln_80.mat'));
il80 = il80.model;

c30 = load(strcat(resultsdir, 'model5/ln_30.mat'));
c30 = c30.model;
c35 = load(strcat(resultsdir, 'model5/ln_35.mat'));
c35 = c35.model;
c40 = load(strcat(resultsdir, 'model5/ln_40.mat'));
c40 = c40.model;
c45 = load(strcat(resultsdir, 'model5/ln_45.mat'));
c45 = c45.model;
c50 = load(strcat(resultsdir, 'model5/ln_50.mat'));
c50 = c50.model;
c55 = load(strcat(resultsdir, 'model5/ln_55.mat'));
c55 = c55.model;
c60 = load(strcat(resultsdir, 'model5/ln_60.mat'));
c60 = c60.model;
c65 = load(strcat(resultsdir, 'model5/ln_65.mat'));
c65 = c65.model;
c70 = load(strcat(resultsdir, 'model5/ln_70.mat'));
c70 = c70.model;
c75 = load(strcat(resultsdir, 'model5/ln_75.mat'));
c75 = c75.model;
c80 = load(strcat(resultsdir, 'model5/ln_80.mat'));
c80 = c80.model;


ages1 = [a30, a35, a40, a45, a50, a55, a60, a65, a70, a75, a80];
ages2 = [i30, i35, i40, i45, i50, i55, i60, i65, i70, i75, i80];
ages3 = [l30, l35, l40, l45, l50, l55, l60, l60, l70, l75, l80];
ages4 = [il30, il35, il40, il45, il50, il55, il60, il65, il70, il75, il80];
ages5 = [c30, c35, c40, c45, c50, c55, c60, c65, c70, c75, c80];
ages_all = cell(1, num_variants);
ages_all{1} = ages1;
ages_all{2} = ages2;
ages_all{3} = ages3;
ages_all{4} = ages4;
ages_all{5} = ages5;

%evolution of selectivity cumulative distribution function with age (combined plots)
for variant = 1:num_variants
    ages = ages_all{variant};
    figure;
    h30 = histogram([ages(1).selectivity.exc; ages(1).selectivity.inh],length(ages(1).selectivity.exc)+length(ages(1).selectivity.inh),'Normalization','cdf');
    h30x = h30.BinEdges(1:(end-1));
    h30y = h30.Values;
    
    h35 = histogram([ages(2).selectivity.exc; ages(2).selectivity.inh],length(ages(2).selectivity.exc)+length(ages(2).selectivity.inh),'Normalization','cdf');
    h35x = h35.BinEdges(1:(end-1));
    h35y = h35.Values;
    
    h40 = histogram([ages(3).selectivity.exc; ages(3).selectivity.inh],length(ages(3).selectivity.exc)+length(ages(3).selectivity.inh),'Normalization','cdf');
    h40x = h40.BinEdges(1:(end-1));
    h40y = h40.Values;
    
    h45 = histogram([ages(4).selectivity.exc; ages(4).selectivity.inh],length(ages(4).selectivity.exc)+length(ages(4).selectivity.inh),'Normalization','cdf');
    h45x = h45.BinEdges(1:(end-1));
    h45y = h45.Values;
    
    h50 = histogram([ages(5).selectivity.exc; ages(5).selectivity.inh],length(ages(5).selectivity.exc)+length(ages(5).selectivity.inh),'Normalization','cdf');
    h50x = h50.BinEdges(1:(end-1));
    h50y = h50.Values;
    
    h55 = histogram([ages(6).selectivity.exc; ages(6).selectivity.inh],length(ages(6).selectivity.exc)+length(ages(6).selectivity.inh),'Normalization','cdf');
    h55x = h55.BinEdges(1:(end-1));
    h55y = h55.Values;
    
    h60 = histogram([ages(7).selectivity.exc; ages(7).selectivity.inh],length(ages(7).selectivity.exc)+length(ages(7).selectivity.inh),'Normalization','cdf');
    h60x = h60.BinEdges(1:(end-1));
    h60y = h60.Values;
    
    h65 = histogram([ages(8).selectivity.exc; ages(8).selectivity.inh],length(ages(8).selectivity.exc)+length(ages(8).selectivity.inh),'Normalization','cdf');
    h65x = h65.BinEdges(1:(end-1));
    h65y = h65.Values;
    
    h70 = histogram([ages(9).selectivity.exc; ages(9).selectivity.inh],length(ages(9).selectivity.exc)+length(ages(9).selectivity.inh),'Normalization','cdf');
    h70x = h70.BinEdges(1:(end-1));
    h70y = h70.Values;
    
    h75 = histogram([ages(10).selectivity.exc; ages(10).selectivity.inh],length(ages(10).selectivity.exc)+length(ages(10).selectivity.inh),'Normalization','cdf');
    h75x = h75.BinEdges(1:(end-1));
    h75y = h75.Values;
    
    h80 = histogram([ages(11).selectivity.exc; ages(11).selectivity.inh],length(ages(11).selectivity.exc)+length(ages(11).selectivity.inh),'Normalization','cdf');
    h80x = h80.BinEdges(1:(end-1));
    h80y = h80.Values;
    close gcf;
    
    figure;
    plot(h30x,100*h30y,'LineWidth',10,'Color',colors(1,:)) % plot young data
    hold on;
    plot(h80x,100*h80y,'LineWidth',10,'Color',colors(11,:)) % plot old data
    
    plot(h35x,100*h35y,'-.','LineWidth',5,'Color',colors(2,:)) % 35
    plot(h40x,100*h40y,'-.','LineWidth',5,'Color',colors(3,:)) % 40
    plot(h45x,100*h45y,'-.','LineWidth',5,'Color',colors(4,:)) % 45
    plot(h50x,100*h50y,'-.','LineWidth',5,'Color',colors(5,:)) % 50
    plot(h55x,100*h55y,'-.','LineWidth',5,'Color',colors(6,:)) % 55
    plot(h60x,100*h60y,'-','LineWidth',5,'Color',colors(7,:)) % 60
    plot(h65x,100*h65y,'-','LineWidth',5,'Color',colors(8,:)) % 65
    plot(h70x,100*h70y,'-','LineWidth',5,'Color',colors(9,:)) % 70
    plot(h75x,100*h75y,'-','LineWidth',5,'Color',colors(10,:)) % 75
    
    %replot the youngest and oldest data to be on top
    plot(h30x,100*h30y,'LineWidth',10,'Color',colors(1,:)) % plot young data
    plot(h80x,100*h80y,'LineWidth',10,'Color',colors(11,:)) % plot old data
    title('Model results','FontSize',30)
    set(gca,'FontSize',25)
    legend('Young network','Old network','Location','Northwest')
    legend('boxoff')
    axis square
    saveas(gcf, strcat(resultsdir, sprintf('cdfs_model%d.fig', variant)))
    export_fig(strcat(resultsdir, sprintf('cdfs_model%d.pdf', variant)), '-painters', '-transparent')
end

%%
%this section is adapted from RunVisionNetworkSimulation.m of E-I Net

%selecting input patches from gratings on which to test network
numOrientations = 4;
numfpa = floor(size(gratings, 3)/numOrientations); %num frames per orientation (angle)
numSamplesPerBatch = a30.numSamplesPerBatch;
gind_r_c = zeros(numSamplesPerBatch, 2); %image indices, row indices, column indices
patchDims = a30.inputDims;
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
        
        %normalize to zero mean and unit variance
        inputData = bsxfun(@minus, inputData, mean(inputData));
        inputData = bsxfun(@times, inputData, 1./std(inputData));
        ginputDataMaster(theta+1,f,:,:) = inputData;
    end
end

%measuring spike counts
spikeCountsExc = zeros(a30.cg_V1e.numCells, numel(ages1));
spikeCountsInh = zeros(a30.cg_V1i.numCells, numel(ages1));
for ag = 1:numel(ages1)
    model = ages1(ag);
    zLength = 50;
    numImages = size(IMAGES, 3);
    numSamplesPerBatch = model.numSamplesPerBatch;
    ind_r_c = zeros(numSamplesPerBatch, 3); %image indices, row indices, column indices
    spikeCountsE = zeros(model.cg_V1e.numCells, 1);
    spikeCountsI = zeros(model.cg_V1i.numCells, 1);
    i_cg_input = FindByName(model.cellGroup, 'input');
    i_cg_output = model.outputCellGroupId;
    for l = 1:numSamplesPerBatch
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
        % normalize to zero mean and unit variance
        inputData = bsxfun(@minus, inputData, mean(inputData));
        inputData = bsxfun(@times, inputData, 1 ./ std(inputData));
        model.snapshot.inputData = inputData;
        if isfield(model, 'finalState')
            % special case simNoReset mode: copy state if after the first iteration
            initialState = model.finalState;
        else
            initialState = struct();
        end
        
        initialState.y{i_cg_input} = inputData * (model.inputScale * model.simTimeStep);
        % simulate the network in batch
        [model, finalState] = NetModel_UpdateFast(model, model.numIterationsPerSample, initialState);
        model.finalState = finalState;
        numSpikesE = sum(finalState.y_history{i_cg_output}, 3);
        numSpikesI = sum(finalState.y_history{i_cg_output+1}, 3);
        for i = 1:model.cg_V1e.numCells
            spikeCountsE(i) = spikeCountsE(i) + numSpikesE(i);
        end
        for i = 1:model.cg_V1i.numCells
            spikeCountsI(i) = spikeCountsI(i) + numSpikesI(i);
        end
    end
    spikeCountsExc(:, ag) = spikeCountsE;
    spikeCountsInh(:, ag) = spikeCountsI;
end
%%
%average spike counts
close gcf;
save('countSpikesExc.txt','spikeCountsExc')
xticks = '30 35 40 45 50 55 60 65 70 75 80';
xticks = split(xticks, ' ');
figure;
distributionPlot(spikeCountsExc, 'xNames', xticks);
xlabel('Age (\# training loops)','FontSize',30,'Interpreter','LaTeX')
ylabel('Relative spike counts','FontSize',30,'Interpreter','LaTeX')
saveas(gcf, strcat(resultsdir, 'excitatory_spike_counts.fig'))
export_fig(strcat(resultsdir, 'excitatory_spike_counts.pdf'), '-painters', '-transparent')

save('countSpikesInh.txt','spikeCountsInh')
figure;
distributionPlot(spikeCountsInh, 'xNames', xticks);
xlabel('Age (\# training loops)','FontSize',30,'Interpreter','LaTeX')
ylabel('Relative spike counts','FontSize',30,'Interpreter','LaTeX')
saveas(gcf, strcat(resultsdir, 'inhibitory_spike_counts.fig'))
export_fig(strcat(resultsdir, 'inhibitory_spike_counts.pdf'), '-painters', '-transparent')


%angles between vectorized input weights
figure;
young = a30;
young_aged_angles = zeros(numel(ages1), young.cg_V1e.numCells);
for i = 1:numel(ages1)
    aged = ages1(i);
    for j = 1:young.cg_V1e.numCells
        young_aged_angles(i, j) = real(acos(dot(young.cellGroup{2}.inputBlock(1).weight(j, :), aged.cellGroup{2}.inputBlock(1).weight(j, :))/(norm(young.cellGroup{2}.inputBlock(1).weight(j, :))*norm(aged.cellGroup{2}.inputBlock(1).weight(j, :)))));
    end
end
young_aged_angles = young_aged_angles';
save('RFangles.txt','young_aged_angles')
distributionPlot(young_aged_angles, 'xNames', xticks);
xlabel('Age (\# training loops)','FontSize',13,'Interpreter','LaTeX')
ylabel('Angle between vectorized young (age 30) and aged input weights','FontSize',13,'Interpreter','LaTeX')
saveas(gcf, strcat(resultsdir, 'angle_btw_vectors.fig'))
export_fig(strcat(resultsdir, 'angle_btw_vectors.pdf'), '-painters', '-transparent')


%lateral and input weight distribution comparisons

% i -> e weights
wei30 = -a30.cellGroup{2}.inputBlock(3).weight;
wei80 = -a80.cellGroup{2}.inputBlock(3).weight;

% i -> i weights
wii30 = -a30.cellGroup{3}.inputBlock(3).weight;
wii80 = -a80.cellGroup{3}.inputBlock(3).weight;

% e -> i weights
wie30 = a30.cellGroup{3}.inputBlock(2).weight;
wie80 = a80.cellGroup{3}.inputBlock(2).weight;

% i -> e weights

figure; histogram(wei30,30,'Normalization','pdf','FaceColor',colors(1,:))
hold on;
histogram(wei80,20,'Normalization','pdf','FaceAlpha',0.5,'FaceColor','b')
xlim([-5 0])
set(gca, 'YScale', 'log','FontSize',25)
xlabel('I$\rightarrow$E weight $W_{EI}$','FontSize',25,'Interpreter','LaTeX')
ylabel('frequency (norm.)','FontSize',25,'Interpreter','LaTeX')
legend('young network', 'old network','Location','Northeast')
legend('boxoff')
ylim([1e-4 10.5])
axis square
saveas(gcf, strcat(resultsdir, 'i_to_e.fig'))
export_fig(strcat(resultsdir, 'i_to_e.pdf'), '-painters', '-transparent')


% i -> i weights

figure; histogram(wii30,20,'Normalization','pdf','FaceColor',colors(1,:))
hold on;
histogram(wii80,20,'Normalization','pdf','FaceAlpha',0.5,'FaceColor','b')
xlim([-3 0])
set(gca, 'YScale', 'log','FontSize',25)
xlabel('I$\rightarrow$I weight $W_{II}$','FontSize',25,'Interpreter','LaTeX')
ylabel('frequency (norm.)','FontSize',25,'Interpreter','LaTeX')
legend('young network', 'old network','Location','Northeast')
legend('boxoff')
ylim([1e-4 10.5])
axis square
saveas(gcf, strcat(resultsdir, 'i_to_i.fig'))
export_fig(strcat(resultsdir, 'i_to_i.pdf'), '-painters', '-transparent')


% e -> i weights
figure; histogram(wie30,32,'Normalization','pdf','FaceColor',colors(1,:))
hold on;
histogram(wie80,20,'Normalization','pdf','FaceAlpha',0.5,'FaceColor','b')
xlim([0 5])
set(gca, 'YScale', 'log','FontSize',25)
xlabel('E$\rightarrow$I weight $W_{IE}$','FontSize',25,'Interpreter','LaTeX')
ylabel('frequency (norm.)','FontSize',25,'Interpreter','LaTeX')
legend('young network', 'old network','Location','Northeast')
legend('boxoff')
ylim([1e-4 10.5])
axis square
saveas(gcf, strcat(resultsdir, 'e_to_i.fig'))
export_fig(strcat(resultsdir, 'e_to_i.pdf'), '-painters', '-transparent')


%threshold distribution comparison
thresh30 = a30.cellGroup{2}.spikeThresh;
thresh80 = a80.cellGroup{2}.spikeThresh;
figure; histogram(thresh80,20,'Normalization','pdf','FaceColor','b','FaceAlpha',0.5)
hold on;
histogram(thresh30,32,'Normalization','pdf','FaceColor',colors(1,:),'EdgeColor','none')
xlim([-0.5 9.5]);
set(gca,'YScale','log','FontSize',25);
xlabel('firing thresholds $\theta$','FontSize',25,'Interpreter','LaTeX')
ylabel('frequency (norm.)','FontSize',25,'Interpreter','LaTeX')
lgt = legend('old network', 'young network','Location','Northeast');
legend('boxoff');
lgt.FontSize = 20;
ylim([1e-4 5.5]);
axis square;
saveas(gcf, strcat(resultsdir, 'thresh.fig'))
export_fig(strcat(resultsdir, 'thresh.pdf'), '-painters', '-transparent')

%input weight distribution comparison and receptive field comparison

Q30 = a30.cellGroup{2}.inputBlock(1).weight;
Q80 = a80.cellGroup{2}.inputBlock(1).weight;

colors = parula(12);
colors(1:12,1:3) = colors([11:-1:1, 12],1:3);

figure;
histogram(Q30(:),32,'Normalization','pdf','FaceColor',colors(1,:))
hold on;
histogram(Q80(:),28,'Normalization','pdf','FaceAlpha',0.5,'FaceColor','b')
set(gca,'YScale','log','FontSize',25)
xlim([-1.1 1.1])
xlabel('input weights $Q$','FontSize',25,'Interpreter','LaTeX')
ylabel('frequency (norm.)','FontSize',25)
legend('young', 'old','Location','Northeast')
legend('boxoff')
ylim([1e-4 10.5])
axis square
saveas(gcf, strcat(resultsdir, 'input_weight.fig'))
export_fig(strcat(resultsdir, 'input_weight.pdf'), '-painters', '-transparent')

figure;
for i=1:4 %picking first four (arbitrary) neurons
    subplot(2,4,i)
    imagesc(reshape(Q30(i,:),8,8)/sqrt(var(Q30(i,:))));
    axis square
    axis off
    colormap('gray')
    
    subplot(2,4,i+4)
    imagesc(reshape(Q80(i,:),8,8)/sqrt(var(Q80(i,:))));
    axis square
    axis off
    colormap('gray')
end

saveas(gcf, strcat(resultsdir, 'young_rfs_top_old_rfs_bottom.fig'))
export_fig(strcat(resultsdir, 'young_rfs_top_old_rfs_bottom.pdf'), '-painters', '-transparent')


%input weight remapping comparison
figure;
r30_35 = testSelectivity(remapWeightsSmooth(a30, a35), ginputDataMaster, numOrientations, numfpa);
r30_40 = testSelectivity(remapWeightsSmooth(a30, a40), ginputDataMaster, numOrientations, numfpa);
r30_45 = testSelectivity(remapWeightsSmooth(a30, a45), ginputDataMaster, numOrientations, numfpa);
r30_50 = testSelectivity(remapWeightsSmooth(a30, a50), ginputDataMaster, numOrientations, numfpa);
r30_55 = testSelectivity(remapWeightsSmooth(a30, a55), ginputDataMaster, numOrientations, numfpa);
r30_60 = testSelectivity(remapWeightsSmooth(a30, a60), ginputDataMaster, numOrientations, numfpa);
r30_65 = testSelectivity(remapWeightsSmooth(a30, a65), ginputDataMaster, numOrientations, numfpa);
r30_70 = testSelectivity(remapWeightsSmooth(a30, a70), ginputDataMaster, numOrientations, numfpa);
r30_75 = testSelectivity(remapWeightsSmooth(a30, a75), ginputDataMaster, numOrientations, numfpa);
r30_80 = testSelectivity(remapWeightsSmooth(a30, a80), ginputDataMaster, numOrientations, numfpa);

r35_30 = testSelectivity(remapWeightsSmooth(a35, a30), ginputDataMaster, numOrientations, numfpa);
r40_30 = testSelectivity(remapWeightsSmooth(a40, a30), ginputDataMaster, numOrientations, numfpa);
r45_30 = testSelectivity(remapWeightsSmooth(a45, a30), ginputDataMaster, numOrientations, numfpa);
r50_30 = testSelectivity(remapWeightsSmooth(a50, a30), ginputDataMaster, numOrientations, numfpa);
r55_30 = testSelectivity(remapWeightsSmooth(a55, a30), ginputDataMaster, numOrientations, numfpa);
r60_30 = testSelectivity(remapWeightsSmooth(a60, a30), ginputDataMaster, numOrientations, numfpa);
r65_30 = testSelectivity(remapWeightsSmooth(a65, a30), ginputDataMaster, numOrientations, numfpa);
r70_30 = testSelectivity(remapWeightsSmooth(a70, a30), ginputDataMaster, numOrientations, numfpa);
r75_30 = testSelectivity(remapWeightsSmooth(a75, a30), ginputDataMaster, numOrientations, numfpa);
r80_30 = testSelectivity(remapWeightsSmooth(a80, a30), ginputDataMaster, numOrientations, numfpa);
close gcf;
figure;
first_remap = [mean(a30.selectivity.exc), mean(r30_35.selectivity.exc), mean(r30_40.selectivity.exc), mean(r30_45.selectivity.exc), mean(r30_50.selectivity.exc), mean(r30_55.selectivity.exc), mean(r30_60.selectivity.exc), mean(r30_65.selectivity.exc), mean(r30_70.selectivity.exc), mean(r30_75.selectivity.exc), mean(r30_80.selectivity.exc)];
second_remap = [mean(a30.selectivity.exc), mean(r35_30.selectivity.exc), mean(r40_30.selectivity.exc), mean(r45_30.selectivity.exc), mean(r50_30.selectivity.exc), mean(r55_30.selectivity.exc), mean(r60_30.selectivity.exc), mean(r65_30.selectivity.exc), mean(r70_30.selectivity.exc), mean(r75_30.selectivity.exc), mean(r80_30.selectivity.exc)];
figure;
plot(loop_number,control,'o','LineWidth',5,'MarkerSize',15,'MarkerFaceColor',colors(11,:),'MarkerEdgeColor','none')
hold on;
plot(loop_number,first_remap,'rs','LineWidth',5,'MarkerSize',15,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none')
plot(loop_number,second_remap,'^','LineWidth',5,'MarkerSize',15,'MarkerFaceColor',colors(6,:),'MarkerEdgeColor','none')
xlabel('Age (\# training loops)','FontSize',30,'Interpreter','LaTeX')
ylabel('Mean selectivity','FontSize',30,'Interpreter','LaTeX')
xlim([30 80]);
ylim([0 1]);
set(gca,'FontSize',25)
lgrf = legend({'Baseline','Young RF + old magnitudes', 'Old RF + young magnitudes'}, 'Location','SouthWest');
legend('boxoff')
set(lgrf,'FontSize',18,'Interpreter','LaTeX');

saveas(gcf, strcat(resultsdir, 'input_weight_remap.fig'))
export_fig(strcat(resultsdir, 'input_weight_remap.pdf'), '-painters', '-transparent')





%young-old parameter swapping

figure;
yng = rmfield(a30, 'response'); %creating a generic copy with response and selectivity fields deleted
yng = rmfield(yng, 'selectivity');
ole = rmfield(a80, 'response');
ole = rmfield(ole, 'selectivity');

oyy = yng;
oyy.cellGroup{2}.inputBlock(1) = a80.cellGroup{2}.inputBlock(1);
oyy.cellGroup{3}.inputBlock(1) = a80.cellGroup{3}.inputBlock(1);
oyy = testSelectivity(oyy, ginputDataMaster, numOrientations, numfpa);
oyy_sel = mean(oyy.selectivity.exc);

yoy = yng;
yoy.cellGroup{2}.inputBlock(3) = a80.cellGroup{2}.inputBlock(3);
yoy.cellGroup{3}.inputBlock(3) = a80.cellGroup{3}.inputBlock(3);
yoy.cellGroup{3}.inputBlock(2) = a80.cellGroup{3}.inputBlock(2);
yoy = testSelectivity(yoy, ginputDataMaster, numOrientations, numfpa);
yoy_sel = mean(yoy.selectivity.exc);

yyo = yng;
yyo.cellGroup{2}.spikeThresh = a80.cellGroup{2}.spikeThresh;
yyo.cellGroup{3}.spikeThresh = a80.cellGroup{3}.spikeThresh;
yyo = testSelectivity(yyo, ginputDataMaster, numOrientations, numfpa);
yyo_sel = mean(yyo.selectivity.exc);

ooy = ole;
ooy.cellGroup{2}.spikeThresh = a30.cellGroup{2}.spikeThresh;
ooy.cellGroup{3}.spikeThresh = a30.cellGroup{3}.spikeThresh;
ooy = testSelectivity(ooy, ginputDataMaster, numOrientations, numfpa);
ooy_sel = mean(ooy.selectivity.exc);

oyo = ole;
oyo.cellGroup{2}.inputBlock(3) = a30.cellGroup{2}.inputBlock(3);
oyo.cellGroup{3}.inputBlock(3) = a30.cellGroup{3}.inputBlock(3);
oyo.cellGroup{3}.inputBlock(2) = a30.cellGroup{3}.inputBlock(2);
oyo = testSelectivity(oyo, ginputDataMaster, numOrientations, numfpa);
oyo_sel = mean(oyo.selectivity.exc);

yoo = ole;
yoo.cellGroup{2}.inputBlock(1) = a30.cellGroup{2}.inputBlock(1);
yoo.cellGroup{3}.inputBlock(1) = a30.cellGroup{3}.inputBlock(1);
yoo = testSelectivity(yoo, ginputDataMaster, numOrientations, numfpa);
yoo_sel = mean(yoo.selectivity.exc);
close gcf;

%bar chart
xlabels = ['YYY'
    'OYY'
    'YOY'
    'YYO'
    'OOY'
    'OYO'
    'YOO'
    'OOO'];

C = [mean(a30.selectivity.exc)
    oyy_sel
    yoy_sel
    yyo_sel
    ooy_sel
    oyo_sel
    yoo_sel
    mean(a80.selectivity.exc)];

figure;
axb = axes;
bp = bar(axb,C,0.5,'FaceColor','flat','LineWidth',5);
ylabel('Mean selectivity','Interpreter','latex','FontSize',20)
ylim([0 1]);
axb.XTickLabel= xlabels;
xtickangle(30);
text(4,0.9,'$Q W \theta =(Y/O,Y/O,Y/O)$','Interpreter','latex','FontSize',20);
axb.FontSize = 20;
set(axb.Legend,'Interpreter','latex','FontSize',20);
bp(1).CData = [1 1 1];
saveas(gcf, strcat(resultsdir, 'young_old_parameter_swaps.fig'))
export_fig(strcat(resultsdir, 'young_old_parameter_swaps.pdf'), '-painters', '-transparent')


%Selectivity when freezing parameter learning rates

input_weight_off = [mean(i30.selectivity.exc), mean(i35.selectivity.exc), mean(i40.selectivity.exc), mean(i45.selectivity.exc), mean(i50.selectivity.exc), mean(i55.selectivity.exc), mean(i60.selectivity.exc), mean(i65.selectivity.exc), mean(i70.selectivity.exc), mean(i75.selectivity.exc), mean(i80.selectivity.exc)];
lateral_weight_off = [mean(l30.selectivity.exc), mean(l35.selectivity.exc), mean(l40.selectivity.exc), mean(l45.selectivity.exc), mean(l50.selectivity.exc), mean(l55.selectivity.exc), mean(l60.selectivity.exc), mean(l65.selectivity.exc), mean(l70.selectivity.exc), mean(l75.selectivity.exc), mean(l80.selectivity.exc)];
both_off = [mean(il30.selectivity.exc), mean(il35.selectivity.exc), mean(il40.selectivity.exc), mean(il45.selectivity.exc), mean(il50.selectivity.exc), mean(il55.selectivity.exc), mean(il60.selectivity.exc), mean(il65.selectivity.exc), mean(il70.selectivity.exc), mean(il75.selectivity.exc), mean(il80.selectivity.exc)];
figure;
plot(loop_number,control,'o','LineWidth',5,'MarkerSize',15,'MarkerFaceColor',colors(11,:),'MarkerEdgeColor','none')
hold on;
plot(loop_number,input_weight_off,'rs','LineWidth',5,'MarkerSize',15,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none')
plot(loop_number,lateral_weight_off,'^','LineWidth',5,'MarkerSize',15,'MarkerFaceColor',colors(4,:),'MarkerEdgeColor','none')
plot(loop_number,both_off,'d','LineWidth',5,'MarkerSize',15,'MarkerFaceColor',colors(7,:),'MarkerEdgeColor','none')

xlabel('Age (\# training loops)','FontSize',30,'Interpreter','LaTeX')
ylabel('Mean selectivity','FontSize',30,'Interpreter','LaTeX')
ylim([0 1]);
xlim([30 80]);
set(gca,'FontSize',25)
lgf = legend({'Baseline','Input weights $Q$ frozen', 'Lateral weights $W$ frozen', 'Both $Q$ \& $W$ frozen'}, 'Location','SouthWest');
legend('boxoff');
set(lgf,'FontSize',18,'Interpreter','LaTeX');

saveas(gcf, strcat(resultsdir, 'freezing_parameter_learning.fig'))
export_fig(strcat(resultsdir, 'freezing_parameter_learning.pdf'), '-painters', '-transparent')


%critical learning period - slowing down learning and extending stimulus
%presentation after maturation

figure;
control = [mean(a30.selectivity.exc), mean(a35.selectivity.exc), mean(a40.selectivity.exc), mean(a45.selectivity.exc), mean(a50.selectivity.exc), mean(a55.selectivity.exc), mean(a60.selectivity.exc), mean(a65.selectivity.exc), mean(a70.selectivity.exc), mean(a75.selectivity.exc), mean(a80.selectivity.exc)];
test = [mean(c30.selectivity.exc), mean(c35.selectivity.exc), mean(c40.selectivity.exc), mean(c45.selectivity.exc), mean(c50.selectivity.exc), mean(c55.selectivity.exc), mean(c60.selectivity.exc), mean(c65.selectivity.exc), mean(c70.selectivity.exc), mean(c75.selectivity.exc), mean(c80.selectivity.exc)];
plot(loop_number,control,'o','LineWidth',5,'MarkerSize',15,'MarkerFaceColor',colors(11,:),'MarkerEdgeColor','none')
hold on;
plot(loop_number,test,'rs','LineWidth',5,'MarkerSize',15,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none')
xlabel('Age (\# training loops)','FontSize',30,'Interpreter','LaTeX')
ylabel('Mean selectivity','FontSize',30,'Interpreter','LaTeX')
ylim([0 1]);
xlim([30 80]);
set(gca,'FontSize',25)
lgf = legend({'Baseline','Slower learning, longer stimulus pres.'}, 'Location','SouthWest');
legend('boxoff');
set(lgf,'FontSize',18,'Interpreter','LaTeX');
saveas(gcf, strcat(resultsdir, 'critical_learning_period_test.fig'))
export_fig(strcat(resultsdir, 'critical_learning_period_test.pdf'), '-painters', '-transparent')

%Gabor character of young vs. old receptive fields as determined by fitting Gabor
%wavelets

%Gabor profiles defined as in Eq. 12 of the paper Zylberberg J, Murphy JT, DeWeese MR (2011). A sparse coding model with synaptically local plasticity and spiking neurons can account for the diverse shapes of V1 simple cell receptive fields. PLoS Computational Biology
%x and y correspond to pixel indices
%Cost function is defined as the squared difference between the gabor profile and the model data, ||G - RF||^2
%fminunc() is used to minimize this cost function and thereby obtain estimates of parameters
young = a30;
old = a80;
young_iw = reshape(young.cellGroup{2}.inputBlock(1).weight', patchDims(1), patchDims(2), []);
old_iw = reshape(old.cellGroup{2}.inputBlock(1).weight', patchDims(1), patchDims(2), []);
young_fits = ones(1, young.cg_V1e.numCells);
old_fits = ones(1, young.cg_V1e.numCells);
x = 1:patchDims(1);
y = 1:patchDims(2);
onevecx = ones(1, patchDims(1));
onevecy = ones(1, patchDims(2));
initvec = [1;0;0;pi;1;1;0.1;0];
gabor_thresh = 0.8; %threshold for cost function value below which fits are deemed good (i.e. receptive field is Gabor-like)
young_unfittable_indices = [];
old_unfittable_indices = [];

for i = 1:young.cg_V1e.numCells
    fprintf('Fitting cell %d out of %d\n', i, young.cg_V1e.numCells);
    RF = young_iw(:, :, i);
    L = @(z) (norm(z(1)*cos(2*pi*z(7)*((x-z(2))'*onevecx*cos(z(4))+onevecy'*(y-z(3))*sin(z(4)))+z(8)).*(exp(-((x-z(2))'*onevecx*cos(z(4))+onevecy'*(y-z(3))*sin(z(4))).^2/2/z(5)^2 - (-(x-z(2))'*onevecx*sin(z(4))+onevecy'*(y-z(3))*cos(z(4))).^2/2/z(6)^2))-RF)^2);
    options = optimoptions('fminunc', 'Display', 'off');
    try
        Gparams = fminunc(L, initvec, options);
        if Gparams(2) < patchDims(1) && Gparams(2) > 0 && Gparams(3) < patchDims(2) && Gparams(3) > 0
            fit_goodness = L(Gparams)/norm(RF)^2;
            if fit_goodness <= gabor_thresh
                young_fits(i) = fit_goodness;
            end
        end
        RF = old_iw(:, :, i);
        L = @(z) (norm(z(1)*cos(2*pi*z(7)*((x-z(2))'*onevecx*cos(z(4))+onevecy'*(y-z(3))*sin(z(4)))+z(8)).*(exp(-((x-z(2))'*onevecx*cos(z(4))+onevecy'*(y-z(3))*sin(z(4))).^2/2/z(5)^2 - (-(x-z(2))'*onevecx*sin(z(4))+onevecy'*(y-z(3))*cos(z(4))).^2/2/z(6)^2))-RF)^2);
        try
            Gparams = fminunc(L, initvec, options);
            if Gparams(2) < patchDims(1) && Gparams(2) > 0 && Gparams(3) < patchDims(2) && Gparams(3) > 0
                fit_goodness = L(Gparams)/norm(RF)^2;
                if fit_goodness <= gabor_thresh
                    old_fits(i) = fit_goodness;
                end
            end
        catch
            old_unfittable_indices = [old_unfittable_indices i];
        end
    catch
        young_unfittable_indices = [young_unfittable_indices i];
    end
end

unfittable_indices = union(young_unfittable_indices, old_unfittable_indices);
total_unfittable = numel(unfittable_indices);
total_fittable = young.cg_V1e.numCells-total_unfittable;
num_gabor_young = numel(find(young_fits ~= 1));
num_gabor_old = numel(find(old_fits ~= 1));

gaborfile = fopen(strcat(resultsdir, 'gabor_fit_results.txt'), 'wt');
line = sprintf('%.2f%% of neurons, in either youth or age, did not have Gabor fits and are excluded from the following statistics', total_unfittable/young.cg_V1e.numCells*100);
fprintf(gaborfile, '%s\n', line);
line = sprintf('%.2f%% of young neurons had good Gabor-like fits, while only %.2f%% of aged neurons did', num_gabor_young/total_fittable*100, num_gabor_old/total_fittable*100);
fprintf(gaborfile, '%s\n', line);
num_staying_gabor = 0; %neurons with receptive fields gabor-like in both youth and old age
num_becoming_gabor = 0; %neurons with receptive fields that are not gabor-like in youth but are in old age
for i = setdiff(1:young.cg_V1e.numCells, unfittable_indices)
    if young_fits(i) ~= 1
        if old_fits(i) ~= 1
            num_staying_gabor = num_staying_gabor + 1;
        end
    else
        if old_fits(i) ~= 1
            num_becoming_gabor = num_becoming_gabor + 1;
        end
    end
end
line = sprintf('%.2f%% of Gabor-like young neurons remained gabor-like in old age, while %.2f%% of neurons that were Gabor-like in youth were rejected as being Gabor-like in old age', num_staying_gabor/total_fittable*100, 100-num_staying_gabor/total_fittable*100);
fprintf(gaborfile, '%s\n', line);
line = sprintf('%.2f%% of neurons were not Gabor-like in youth but became Gabor-like in old age', num_becoming_gabor/total_fittable*100);
fprintf(gaborfile, '%s\n', line);
fclose(gaborfile);
type(strcat(resultsdir, 'gabor_fit_results.txt'))

%Dependence between receptive field similarity and strength of mutual
%inhibition of pairs of excitatory neurons
model = a30;
ee = model.cellGroup{2}.inputBlock(3).weight*model.cellGroup{3}.inputBlock(2).weight;
overlap = zeros(model.cg_V1e.numCells*(model.cg_V1e.numCells-1)/2,2);
inhibition = zeros(model.cg_V1e.numCells*(model.cg_V1e.numCells-1)/2,2);
c=0;
for i = 1:model.cg_V1e.numCells
    for j = (i+1):model.cg_V1e.numCells
        c=c+1;
        similarity = dot(model.cellGroup{2}.inputBlock(1).weight(i, :), model.cellGroup{2}.inputBlock(1).weight(j, :))/(norm(model.cellGroup{2}.inputBlock(1).weight(i, :))*norm(model.cellGroup{2}.inputBlock(1).weight(j, :)));
        overlap(c, 1) = similarity;
        overlap(c, 2) = similarity;
        inhibition(c, 1) = ee(i, j);
        inhibition(c, 2) = ee(j, i);
    end
end

figure;
scatter(inhibition(:), overlap(:)) 
hold on; plot(0:40,0*(0:40),'k--','linewidth',5)
xlabel('$|W_{EI} W_{IE}|$','Interpreter','latex','fontsize',20)
ylabel('RF overlap','Interpreter','latex','fontsize',20)
set(gca,'fontsize',20)
% saveas(gcf, strcat(resultsdir, 'rf_similarity_vs_inhibition.fig'))
% export_fig(strcat(resultsdir, 'rf_similarity_vs_inhibition.pdf'), '-painters', '-transparent')
end


