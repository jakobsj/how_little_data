%
% Example script of running a parameter sweep using parobj
%
% Jakob S. Joergensen (jakj@dtu.dk), 2014.

% Initialize an empty parameter object po:
po = parobj

% Make three parameters, of different types and number of values
po.setNames({'height','weight','age'})
po.setTypes({'%1.2f','%2.1f','%d'})
po.setValues({[1.69, 1.84, 1.79], [67, 82, 75], [21, 40, 59]})

% Form the array of all possible parameter sets in po.array
po.buildArray()

% The parameter object is now ready for sweeping since the array of all
% parameter sets is created, which can be seen by displaying po:
po

% Sweep over all sets of parameters for the simfunc. Second input is
% boolean: true if the simfunc accepts one extra input argument compared to
% the number of parameters (here, 3, so a total of 4 inputs arguments),
% where the 4th will be a filename automatically generated from the
% parameter values for saving files systematically. false if the simfunc
% only accepts the number of sweeping parameters as input (here 3).
simfunc = 'printdata';
do_pass_filename = false;
po.parsweep(simfunc,do_pass_filename)