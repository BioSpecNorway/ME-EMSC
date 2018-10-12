function [correctedSpectraForIteration, residualsFromIteration, ...
    parameters, numberOfIterations] = iterationSteps(specrumRowForCorrection, spectraNumber, correctedSpectraForIteration, ...
                                                        EMSCScaleModel, maxIterationNumber, weights, wn, alpha0, gamma, ...
                                                        options, PCnumber )

    RMSE = zeros(1, options.maxIterationNumber);
    for iterationNumber = 2 : maxIterationNumber   
        % Scale the reference spectrum in each iteration
        if options.scaleRef 
            [correctedSpectraForIteration, ~, ~] = basicEMSCsolver(correctedSpectraForIteration, EMSCScaleModel); 
        end 

        % Input spectrum for calculating the imaginary part of the refractive index, weights applied
        correctedSpectraForIteration = correctedSpectraForIteration .* weights; 
        % Reference spectrum for EMSC, weights applied 
        referenceSpectrum = correctedSpectraForIteration; 

        % Input spectrum for calculating the imaginary part of the refractive index, negative parts set to zero
        correctedSpectraForIteration( correctedSpectraForIteration < 0 ) = 0; 
        
        % Reference spectrum for EMSC, negative parts set to zero
        if options.PositiveRefSpectrum  
            referenceSpectrum = correctedSpectraForIteration; 
        end  
        
        % create ME-EMSC model
        [MieEMSCmodel] = make_ME_EMSCmodel(referenceSpectrum, correctedSpectraForIteration, wn, alpha0, gamma, options, PCnumber); 
        % Calculate corrected spectrum and parameters from EMSC
        [correctedSpectraForIteration, residualsFromIteration, parameters] = ME_EMSCsolver(specrumRowForCorrection, MieEMSCmodel); 
        
        %% Calculate root mean square error 
        if options.fixIterationNumber == false 
            RMSE(iterationNumber) = round( sqrt( ( 1 / ( size( specrumRowForCorrection, 2 ) ) ) * sum( ( residualsFromIteration ) .^ 2 ) ), 4);
        end
        %% Stop criterion
        if iterationNumber == options.maxIterationNumber
            msg = ['Spectrum no. ' num2str(spectraNumber) ': Number of iterations (maxIterationNumber) should be bigger.'];
            warning(msg)
            numberOfIterations = options.maxIterationNumber;
        elseif options.fixIterationNumber && iterationNumber < options.fixIterationNumber
            continue;
        elseif iterationNumber == options.fixIterationNumber
            numberOfIterations = iterationNumber;
            break; 
        elseif iterationNumber > 2 && options.fixIterationNumber == false 
            if RMSE(iterationNumber) == RMSE(iterationNumber-1) && RMSE(iterationNumber) == RMSE(iterationNumber-2) || RMSE(iterationNumber) > RMSE(iterationNumber-1) % When RMSE has been stable up to four digits for two iterations, or it increases   
                numberOfIterations = iterationNumber;  
                break;
            end           
         end 
    end
end

