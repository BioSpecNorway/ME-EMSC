function[correctedSpectra, residuals, EMSCparameters, numberOfIterations, options] = ME_EMSC(normalizedReferenceSpectrum, spectraForCorrection, wn, options)
    %  Correcting Mie scattering in infrared spectra.  
    %
    %  ---------------------------------------------------------------------------------
    %  Written by: 
    %  Johanne Solheim, Evgeniy Gunko, Dennis Peterson, Tatiana Konevskikh, Achim Kohler                                             
    %                                                                                               
    %  Faculty of Science and Technology (REALTEK)                                   
    %  Norwegian Unversity of Life Sciences (www.nmbu.no)                                                                                                                             
    %                                                                                               
    %  Post address:                                                                                                                                                                         
    %  PO Box 5003, 1432 Aas, Norway                                                               
    %                                                                                               
    %  First version: 2016                                                                      
    %  Revision: 2018  
    %  ---------------------------------------------------------------------------------
    %
    %  Input: 
    %  normalizedReferenceSpectrum  - Reference spectrum used in 1st iteration, row vector
    %  pureAbsorbanceSpectra        - Updated 'pure absorbance' spectrum, row vector 
    %  specraForCorrection          - specrumRowForCorrection spectra to correct, one spectrum in each row
    %  wn                           - Wavenumbers corresponding to all spectra, coloumn vector 
    %  options                      - Parameters for correcton.
    %  weightFunctionParams         - parameters for the weight function
    %  
    %  Output: 
    %  correctedSpectra   - Corrected spectra 
    %  residuals          - Residuals for every correction, each row corresponds to one spectrum  
    %  EMSCparameters     - EMSC parameters for every correction, each row corresponds to one spectrum. 
    %                       EMSC parameters are given in the order described in the ME_EMSCsolver function
    %  numberOfIterations - Number of iteration for each spectrum, row vector 

    %% Initialization with default options
    % Default options for the corrections 
    options_default.maxIterationNumber = 45; % Int: Maximal number of iterations 
    options_default.scaleRef = true;
    options_default.PCnumber = false; % Int/False: Number of loadings used in the EMSC-model. Overwrites ExplainedVariance if both are set. 
    options_default.PositiveRefSpectrum = true; % True/False: Set negative parts of the reference spectrum to zero 
    options_default.fixIterationNumber = false; % Int: Fixed number of iterations 
    options_default.mode = 'Correction'; 
    options_default.Weights = true; 
    
    options_default.Weights_InflectionPoints = {[3700 2550], [1900 0]}; % Inflection points of weight function, decreasing order 
    options_default.Weights_Kappa = {[1 1], [1 0]}; % Slope of weight function at corresponding inflection points (options.Weights_InflectionPoints)

    options_default.minRadius = 2;
    options_default.maxRadius = 7.1;
    options_default.minRefractiveIndex = 1.1;
    options_default.maxRefractiveIndex = 1.4;
    options_default.samplingSteps = 10;

    options_default.radius = linspace(options_default.minRadius, options_default.maxRadius, options_default.samplingSteps); % radius of the cell
    options_default.n_zero = linspace(options_default.minRefractiveIndex, options_default.maxRefractiveIndex, options_default.samplingSteps); % constant part of the refractive index
    options_default.h = 0.25; % scale index in gamma range
    options_default.ExplainedVariance = 99.96; % Float/False: Explained varance in the set of Mie extinction curves. Used to determine number of loadings in the EMSC-model.  
    
    options_default.plotResults = false; 
    
    %% Overwrite settings 
    % Overwrite default settings if specified 
    options_all=fieldnames(options_default);
    for n=1:length(options_all)
        if ~isfield(options,options_all{n})
            options = setfield(options,options_all{n},getfield(options_default,options_all{n}));
        end
    end
    
    options.radius = linspace(options.minRadius, options.maxRadius, options.samplingSteps); % radius of the cell
    options.n_zero = linspace(options.minRefractiveIndex, options.maxRefractiveIndex, options.samplingSteps); % constant part of the refractive index
   
    
    % Overwrite maxIterationNumber if fixIterationNumber is sepcified       
    maxIteraionNumber = options.maxIterationNumber;  
    if options.fixIterationNumber
        maxIteraionNumber = options.fixIterationNumber;
    end
    
    %% Errors and warnings
    % Raise error message if mode is incorrect
    if not(strcmp(options.mode, 'PreRun')) && not(strcmp(options.mode, 'Correction'))
        msg = ['Choose options.mode either ' char(39) 'Correction' char(39) ' or ' char(39) 'PreRun' char(39)];
        error(msg)  
    end 
    
    % Raise warning if only one iteration is chosen 
    if options.maxIterationNumber == 1
        msg = 'Using only one iteration will not result in a proper correction.';
        warning(msg)
    end 
        
    %% Settings for Parameter optimization mode 
    if strcmp(options.mode, 'PreRun')
        options.PositiveRefSpectrum = false; 
        options.Weights = false; 
        options.plotResults = true; 
        if size(spectraForCorrection,1)>40
            msg = ['Number of spectra exceeds 40 in options.mode ' char(39) 'PreRun' char(39) '. ' ...
                newline 'A fewer number of spectra is adviced for this mode. Do you wish to proceed? [Y/N]: '];
            answer = input(msg, 's'); 
            if strcmp(answer, 'N')
                return 
            end 
        end 
    end 
    
    %% Build weights for the reference spectrum 
    if options.Weights 
        weights = calcWeightFunction(wn, options);
    else 
        weights = ones(1,length(wn));
    end 
    
    %% Calculate alpha0 and gamma 
    % physical parameters for calculation of Mie extinction curves
    alpha0 = ( 4 * pi * options.radius .* ( options.n_zero - 1 ) ) * (1e-6); 
    gamma = options.h * log(10) ./ ( 4 * pi * 0.5 * pi * ( options.n_zero - 1 ) .* options.radius * (1e-6));
    
    %% Initialize the correction
    numberOfSpectra = size(spectraForCorrection,1); % number of spectra to be corrected 
    
    numberOfIterations = ones(1, numberOfSpectra);
    
    %% first iteration
    % Set pureAbsorbanceSpectra equal to normalizedReferenceSpectrum
    pureAbsorbanceSpectrum = normalizedReferenceSpectrum; 
     
    % create ME EMSC model
    [MieEMSCmodelForFirstIteration, ...
        PCnumber] = make_ME_EMSCmodel(normalizedReferenceSpectrum, pureAbsorbanceSpectrum, wn, alpha0, gamma, options);
    
    % EMSC model for scaling the reference spectrum in each iteration
    if options.scaleRef && maxIteraionNumber > 1
        [EMSCScaleModel] = make_basicEMSCmodel(normalizedReferenceSpectrum,wn,1); 
    else
        EMSCScaleModel = 0;
    end 
    
    % Correction and EMSC parameters from the first iteration
    [correctedSpectra, ...
        residuals, ...
        EMSCparameters] = ME_EMSCsolver(spectraForCorrection, MieEMSCmodelForFirstIteration); 
    
    %% iterations loop  
    if maxIteraionNumber > 1
        for spectrumNumber = 1 : numberOfSpectra  
            [correctedSpectraForIteration, residualsFromIteration, ...
                parameters, numberOfIterationsForSpectra ] = iterationSteps( spectraForCorrection(spectrumNumber, :), spectrumNumber, ...
                                                                correctedSpectra(spectrumNumber, :), EMSCScaleModel, ...
                                                                maxIteraionNumber, weights, wn, alpha0, gamma, ...
                                                                options, PCnumber );
            correctedSpectra(spectrumNumber,:) = correctedSpectraForIteration;
            EMSCparameters(spectrumNumber,:)=parameters;
            residuals(spectrumNumber,:) = residualsFromIteration;    
            numberOfIterations(1, spectrumNumber) = numberOfIterationsForSpectra;
        end        
    end