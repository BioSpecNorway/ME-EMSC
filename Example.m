%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                               %
%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                 %
%            %                                %                                                 %
%            %        --- EXAMPLE ---         %                                                 %
%            %                                %                                                 %
%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                 %
%                                                                                               %
%                                                                                               %
%  Johanne Solheim, Evgeniy Gunko, Achim Kohler                                                 %
%                                                                                               %
%  Faculty of Science and Technology (REALTEK)                                                  %  
%  Norwegian Unversity of Life Sciences (www.nmbu.no)                                           %
%                                                                                               %
%  Post address:                                                                                %
%                                                                                               %
%  PO Box 5003, 1432 Aas, Norway                                                                %
%                                                                                               %  
%                                                                                               %
%  Description: Example of how to change the input parameters for the ME-EMSC correction.       %
%               Parameters are changed in the third and fourth block.                           %
%                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('data'))
addpath(genpath('computing'))
addpath(genpath('helpers'))
addpath(genpath('plots'))

clear all;
close all;

%% Load the data set from the path
load 'data/MatrigelSpectrum.mat'; % Matrigel spectrum (Reference)
load 'data/measuredSpectra.mat'; % Measured spectrum for correction 

%% Set options 
% Choose mode 'Correction' for correcting the whole data set, 
% and mode 'PreRun' for optimizing the input parameters 
options.mode = 'Correction';      

% Change parameter A, number of loadings (principal component directions)
% included in the ME-EMSC model. This is done in two different ways: 
% 1) Set number of loadings directly by specifying options.PCnumber (default false). 
% By specifying this parameter, options.ExplainedVariance is overwritten 
options.PCnumber = 12; 
% 2) Change the limit of the explained variance captured by the loadings (default 99.96). 
% NB! When options.PCnumber is specified, this limit is not used. 
options.ExplainedVariance = 99.99; 

% Change the weight function: 
% To turn of weighting, set options.Weights = false; (default true) 
% To change the appearance of the weight function you can change inflection points. 
% Inflection points are given in decreasing order (default {[3700 2550], [1900 0]}): 
options.Weights_InflectionPoints = {[3800 2700], [1800 0]}; 
% (The last element can be changed from 0 (no downweighting in the end of the spectrum) 
% if downweighting in the end is desired. Example is not shown). 
% Slope of tangent hyperbolic functions at corresponding inflection points are changed 
% by specifying options.Weights_Kappa (default {[1 1], [1 0]}): 
options.Weights_Kappa = {[0.5 0.4], [1 0]};
% By turning on the plotting, you see the effect of cahnging the weight function.  

% Turn on plotting of results: 
% For mode 'Correction' the plotting is default set to false (no plotting).
% For mode 'PreRun' the plotting is default set to true (plot results). 
% Turn on plotting by changing options.plotResults: 
options.plotResults = true; 

% Change maximum number of iterations (default 45): 
% If maximum number of iteration is reached for a spectrum, the correction terminates, 
% and a warning is raised. 
options.maxIterationNumber = 15;

% Set a fixed number of iterations (defualt false): 
% All spectra will be corrected for the specified fixed number of iterations. 
% This option overwrites the regular stop criterion and options.maxIterationNumber.  
% Correct all spectra with 3 iterations by setting: options.fixIterationNumber = 3;

% Change physical parameters: 
% Physical parameters (radius and constant part of the real part of the refractive index) 
% are specified by setting the upper and lower limit of the range. Ranges consist of 10 equidistant points. 
options.minRadius = 3; % Default: 2 
options.maxRadius = 9.1; % Default: 7.1 
options.minRefractiveIndex = 1.2; % Default: 1.1 
options.maxRefractiveIndex = 1.5; % Default: 1.4  

%% Simple quality test - determines which spectra to discard based on a limit for RMSE
RMSE_limit = 1000; % Upper limit for RMSE. Default: infinity (1000)     

%% Data converting
referenceSpectrum = Mat_1000_4000(:,2)'; % Reference spectrum (Matrigel), row vector 
normalizedReferenceSpectrum = referenceSpectrum/max(referenceSpectrum); % Normalize reference spectrum 
wn_ref = Mat_1000_4000(:,1); % Wavenumbers, coloumn vector 
measuredSpectra = Spectra.d; % Spectra to be corrected, one spectrum per row
wn_raw = str2num(Spectra.v); % Wavenumbers, coloumn vector 

%% Adjust the wavenumber region and values of the reference spectrum to be compatible with the raw dataset
[normalizedReferenceSpectrum, measuredSpectra, wn] = adjustWavenumbers(normalizedReferenceSpectrum, wn_ref, measuredSpectra, wn_raw); 

%% Selected spectra for correction 
selectedSpectraNumbersForCorrection  = [1:size(measuredSpectra,1)]; % Select all spectra 
selectedSpectraForCorrection = measuredSpectra(selectedSpectraNumbersForCorrection, :); % Selection of spectra for correction 

%% Run Mie correction 
[correctedSpectra, residuals, EMSCparameters, numberOfIterations, options] = ME_EMSC(normalizedReferenceSpectrum, selectedSpectraForCorrection, wn, options);

%% Calculate RMSE for all spectra for quality test 
for i=1:length(selectedSpectraNumbersForCorrection)
    RMSE(i) = sqrt((1/(size(selectedSpectraForCorrection,2)))*sum((residuals(i, :)).^2)); 
end 


%% Remove spectra with RMSE > RMSE_limit 
ProcessedQT = nan(size(correctedSpectra));
ResidualsQT = nan(size(residuals));
EMSCParametersQT = nan(size(EMSCparameters)); 
numberOfIterationsQT = nan(size(numberOfIterations)); 
ProcessedQT(RMSE<RMSE_limit,:) = correctedSpectra(RMSE<RMSE_limit,:); 
ResidualsQT(RMSE<RMSE_limit,:) = residuals(RMSE<RMSE_limit,:);
EMSCParametersQT(RMSE<RMSE_limit,:) = EMSCparameters(RMSE<RMSE_limit,:);
numberOfIterationsQT(RMSE<RMSE_limit) = numberOfIterations(RMSE<RMSE_limit); 
discardedSpectraNumber = find(RMSE>RMSE_limit); 

%% Plot results 
if options.plotResults
    % Plot of reference spectrum with weights
    referenceSpectrumPlot(normalizedReferenceSpectrum, wn, options);

    % Plot the corrected spectra from the last iteration step 
    correctedSpectrumPlot(ProcessedQT, normalizedReferenceSpectrum, wn)

    % Plot RMSE values and RMSE limit
    RMSEvaluesPlot(RMSE, RMSE_limit, selectedSpectraNumbersForCorrection);

    % Plot residuals at the last iteration step for every corrected spectrum
    residualsPlot(ResidualsQT, wn)
end 









