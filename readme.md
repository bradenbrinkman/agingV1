Reference paper: Talyansky & Brinkman, "Dysregulation of excitatory
neural firing replicates physiological and functional changes in aging
visual cortex," *PLoS Computational Biology* (2020)

***Novel scripts developed for this project***

The MATLAB scripts in this package were written for MATLAB version
2020b, but they should also work with at least MATLAB versions 2019b and
2018b. For full use the Curve Fitting and Optimization Toolboxes must be
installed. The folders “cm1,” “cm7,” and “cm16” in this package contain
the results for movies 01, 07, and 16 from the [*CatCam
database*](https://zenodo.org/record/46481#.X8pfbKpKiu4). 
Download the corresponding movie from this database. (The processing 
of these movies is done by the provided mainAnalyzer.m script). 
Each of the cm folders contains the results for 5 variations of the model 
(model1, model2, etc.), indexed as follows:

1.  Normal aging

2.  Input weights frozen after youth

3.  Lateral weights frozen after youth

4.  All weights frozen after youth

5.  Critical learning (after youth, rate reduced by 10 and stimulus
     presentation time extended by 10)

**mainAnalyzer.m**

-   Driver function that reproduces the results reported in the paper

    -   comment out the line that calls NetworkAger.m if you just want 
        to reproduce the figures in the paper, or else it will retrain
        the entire network and replace the model data files provided
        for the cat movie inputted

    -   This script will not produce most of the Supplementary figures,
        as these involve running the network for multiple different
        initial training movies or different initial parameters.

        -   Supplementary Figure 10 is the only supplementary figure
            reproduced by this script (RF similarity versus |W\_EI W\_IE|)

        -   To reproduce results for different target spike rate p\_E(0)
            (Supplementary Figure 11), please modify line 74 of
            networkAger.m

        -   To reproduce results for different initial movies
            (Supplementary Figures 12-19), please modify the movie
            argument (“naturalmovie”) in your call of mainAnalyzer.m

-   Has sections adapted from code from Olshausen & Field (1997) and
    from the function RunVisionNetworkSimulation.m in the original E-I
    Net package (see “Dependencies”). The original scripts are
    excluded from this package because their functionality is
    encompassed by mainAnalyzer.m and networkAger.m

-   Requires Curve Fitting Toolbox for successful call of
    remapWeightsSmooth.m and Optimization Toolbox for Gabor-ness
    analysis

**networkAger.m**

-   function for training E-I Net

-   Trains five variations of the model in parallel (described above)

-   Has sections adapted from the functions demo3\_TrainEINetFast.m and
    RunVisionNetworkSimulation.m in the original E-I Net package (see
    “Dependencies”). These functions are excluded from this package
    because their functionality is encompassed by mainAnalyzer.m and
    networkAger.m and therefore they are not used.

**orientationSelAnalysis.m**

-   Helper function for computing orientation selectivity from neuron
    responses to the grating stimuli

-   Uses the file gratings.mat, containing moving grating stimuli
    generated with VisualStimulusToolbox (see “Dependencies”), a
    procedure also included in mainAnalyzer.m gratings.mat.

**remapWeightsSmooth.m**

-   Helper function for computing remapped connection weight matrices by
    remapping described in Equation (12) of the paper

-   Requires Curve Fitting Toolbox

**testSelectivity.m**

-   Helper function for testing the selectivity of a network

-   Adapted from E-I Net’s RunVisionNetworkSimulation.m (see
    “Dependencies”). This original script is excluded from this
    package because its functionality is encompassed by other
    functions

***Dependencies and base code upon which this package is built***

This project’s code builds on the E-I Net model referenced above and
linked to below, in addition to using code from other external packages
for presenting grating stimuli and making figures. The necessary code
from these packages is included here for user convenience, and is
subject to the licenses of the original packages.

All the files in the neurosim and util subfolders are preserved from the
E-I Net package, while distributionPlot, export\_fig\_container, and
VisualStimulusToolbox-master contain the necessary code and/or license
files from these packages.

**Original E-I Net package**

-   Paper: King et al., [*“Inhibitory interneurons decorrelate
    excitatory cells to drive sparse code formation in a spiking model
    of V1,” *Journal of Neuroscience** *33 (13) 5475-5485 (2013)*](https://www.jneurosci.org/content/33/13/5475.long).

-   Code: [*https://github.com/paulking/ei\_net*](https://github.com/paulking/ei_net)

-   Note: some component files of the E-I Net package not essential to
    this project have been removed from this package

-   See also SAILnet, the original model upon which E-I Net was
    developed.

    -   Paper: [*“A Sparse Coding Model with Synaptically Local
        Plasticity and Spiking Neurons Can Account for the Diverse
        Shapes of V1 Simple Cell Receptive Fields,” PLoS Computational Biology, 7(10): e1002250 (2011).*](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002250)

    -   Code:
        [*http://jzlab.org/sailcodes.html*](http://jzlab.org/sailcodes.html)

**Supplementary external packages**

-   distributionPlot:  [*https://www.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m*](https://www.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m)

    -   Note: the violin plots provided in Talyansky & Brinkman (2020)
        are generated with the Jupyter notebook
        ViolinPlots\_TalyanskyBrinkman.ipynb included in the main
        directory (requires Python 3.6+), which uses the data files
        countSpikesExc.txt, countSpikesInh.txt, and RFangles.txt. The
        package distributionPlot is included only for the convenience
        of generating similar plots directly in MATLAB.

-   export\_fig: [*https://www.github.com/altmany/export_fig*](https://www.github.com/altmany/export_fig)

-   VisualStimulusToolbox: [*https://zenodo.org/record/154061\#.X8FP0hNKh0s*](https://zenodo.org/record/154061#.X8FP0hNKh0s)


