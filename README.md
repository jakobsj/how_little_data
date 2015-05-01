# how_little_data
This code accompanies the journal article:

How little data is enough? Phase-diagram analysis of sparsity-regularized X-ray computed tomography
Jakob S. Joergensen and Emil Y. Sidky
Phil. Trans. R. Soc. A, Vol. 373, 20140387.
http://dx.doi.org/10.1098/rsta.2014.0387

# Contact:
Jakob S. Joergensen, jakj@dtu.dk

Emil Y. Sidky, sidky@uchicago.edu

The package contains code to run all simulations presented in the article and process all data into resulting data files. 

# Data and code to make figures
The resulting data used for the article as well as scripts to produce the figures of the paper are available in Dryad Digital Repository at
http://dx.doi.org/10.5061/dryad.3jg57

#Contributors

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


# Structure of code
The code is structured in the following way:

The directory "phasediagrams_code" contains code to run all phase diagram simulations and process results into phase diagram results ready for display. The directory "largescale" contains code and data to run all large-scale simulations including the walnuts reconstructions. See READMEs in each directory for further details.

