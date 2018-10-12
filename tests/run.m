%% This script is used to run test functions

addpath(genpath('../data'))
addpath(genpath('../computing'))
addpath(genpath('../helpers'))
addpath(genpath('../tests'))

clear all;
close all;

mie_hulst_complex_scaledTest();
mieCorrectionTest();
disp('Tested'); 