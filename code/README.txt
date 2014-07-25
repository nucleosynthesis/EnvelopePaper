##################################################
#		Instructions for running bias studies        #
##################################################

For recommended use skip to part 2) Running on lxbatch

1) Running locally

The toy throwing and fitting is done in c++ through RooFit
with two classes PdfModelBuilder.cc and ProfileMultiplePdfs.cc
These live in the src directory with their headers in interface.

First you need to setup the environment:

[shell:] source setup.sh # should be name changed to this
[shell:] source setup_root_lxplus.sh

Then compile:

[shell:] make

The bias studies get run by a script called:

bias_study.py

This requires a config file like the one in:

dat/bias_study_config.dat

The config gets loaded by the class:

python/config.py

The output of this is saved in a root file as a TGraph which is the profile likelihood scan
of the fit of each model to each toy. The absolute likelihood is not saved on the delta from
the best fit. If the option savePVal=True then it will also save another TGraph with extension
"_correction". This saves a conversion of the absolute likelihood to an equivalent chi2 with
no free degrees of freedom.

You can then compute the bias using the class python/biasComputer.py

It is CPU intensive to run lots of models with lots of toys so this can be automated via the
lxbatch system. This is the recommeneded way of running.

2) Running on lxbatch

The whole set of jobs should be configured through a datfile similar to:

dat/config_1storders.dat
dat/config_allorders.dat

You can submit and then post-process jobs on the batch system using some stuff in the scripts
directory.

2a) Submission (this jobs take a while so need a long queue)

./scripts/sub_bias.py -d dat/config_1storders.dat -q 8nh

2b) Post processing - bias extraction (this jobs take a short time so need a short queue)

./scripts/extract_bias.py -d dat/config_1storders.dat -s -q 8nm

2c) Making the summary file

./scripts/make_bias_summary_file.py -f BiasSummary.root -d dat/config_1storders.dat

The output file should now contain graphs for the pull mean, pull width and the coverage.

These can be plotted using the scripts:

plotbiasstudies.py
plotcoverages.py
