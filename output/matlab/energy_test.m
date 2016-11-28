clear all;
close all;

epsilon=load('../refinement_epsilon.dat');
derivative_error=load('../refinement_derivative_error.dat');
force_error=load('../refinement_force_error.dat');
plot(log(epsilon),log(derivative_error))
figure
plot(log(epsilon),log(force_error))