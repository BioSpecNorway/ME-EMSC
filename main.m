%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                               %
%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                 %
%            %                                %                                                 %
%            %           --- MAIN ---         %                                                 %
%            %                                %                                                 %
%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                 %
%                                                                                               %
%                                                                                               %
%  Johanne Solheim, Evgeniy Gunko, Dennis Petersen, Tatiana Konevskikh, Achim Kohler            %
%                                                                                               %
%  Faculty of Science and Technology (REALTEK)                                                  %  
%  Norwegian Unversity of Life Sciences (www.nmbu.no)                                           %
%                                                                                               %
%  Post address:                                                                                %
%                                                                                               %
%  PO Box 5003, 1432 Aas, Norway                                                                %
%                                                                                               %
%  First version: 2016                                                                          %
%                                                                                               %
%  Revision: 2018                                                                               %
%                                                                                               %
%  Description: Mie Extinction EMSC (ME-EMSC) correction for a data set                         %
%                                                                                               %                                                                                              %                                                                             %
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
options.mode = 'PreRun';  % Choose between: 'PreRun' and 'Correction'   
% Mode 'PreRun' is used for optimizing the input parameters. Use only a selection of 10-20 spectra for the pre run. 
% Mode 'Correction' is used for running the Mie correction on the whole
% data set.
% See Example.m for how to change the default options.

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
% selectedSpectraNumbersForCorrection  = [1:size(measuredSpectra,1)]; % All spectra are used 
selectedSpectraNumbersForCorrection = [1:5]; % Only use a selection of the spectra 
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