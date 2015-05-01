function [] = printdata(h,w,a)
%PRINTDATA Example function for running parameter sweep. Print data.
%
% Jakob S. Jorgensen (jakj@dtu.dk), 2014.

fprintf('Height : %.2f m\n',h)
fprintf('Weight : %.1f kg\n',w)
fprintf('Age    : %d years\n',a)
fprintf('BMI    : %.2f\n', w/h^2)