Software needed to generate the phase diagram results in:

How little data is enough? Phase-diagram analysis of sparsity-regularized X-ray computed tomography
Jakob S. Joergensen and Emil Y. Sidky
Phil. Trans. R. Soc. A, Vol. 373, 20140387.
http://dx.doi.org/10.1098/rsta.2014.0387

Contributors to code used in the study (besides the authors):

AIR Tools
Maria Saxild-Hansen, Per Christian Hansen, Jakob S. Jørgensen.
Version 1.0.
http://www2.compute.dtu.dk/~pcha/AIRtools/

export_fig
Yair Altman, Oliver J. Woodford
http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig

SNOWMAKER
Michael B. McCoy
http://www.maths.manchester.ac.uk/~mlotz/SNOWMAKER.tar.gz

Tabulated polytope phase-transition data
Jared Tanner
http://people.maths.ox.ac.uk/tanner/data/polytope.mat

get_test_image_isotv_qr_repeat
Christian Kruschel, Jakob S. Jørgensen

mosek_wrap
Anders Skajaa, Jakob S. Jørgensen


The code is structured in the following way:

- code/exec contains the main scripts to execute simulations for Donoho-Tanner (dt) and Amelunxen-Lotz-McCoy-Tropp (almt) phase diagrams as well as process the data into resulting phase diagrams ready for display.
- code/ext contains the various external packages needed.
- code/load_funcs contains functions for loading in simulation result files to be used in the processing into phase diagram data.
- code/run_funcs contains wrapper function to simplify the call to the simulation code for the various image class and measurement setups. 
- code/sim_funcs contains the core simulation functions as well as functions needed in the processing of simulation results into phase diagram data.
- data_raw/ is where the simulation results will be saved when running the execution scripts.
- data_processed_user/ is where the final processed phase diagram files will be deposited when running the processing scripts, which require a complete set of simulation results to be present in the data_raw directory.

Please note that most of the simulations take a very long time to complete, for example each phase diagram requires solving on the order of 100000 individual optimization problems. The supplied code executes the simulations on a single machine, which will most likely not be feasible for completing the full simulation. The full set of simulations for the article were run on a cluster; code to do this is not included. One may change the parameters in the execution and processing scripts to produce smaller data sets faster.

Questions or comments:
Send message to Jakob S. Joergensen: jakj@dtu.dk

