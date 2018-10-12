function mieCorrectionTest()
    %% Test for checking correction with known data
    
	load 'data/MieCorrection/input_dataset_1.mat';
	load 'data/MieCorrection/output_dataset_1.mat';

    %% convert input date into new format
    ZProcessedOriginal = zeros(60, 1556);
    for i = 1:60
        ZProcessedOriginal(i, :) = ZProcessed(i).d;
    end

    ZResidualsOriginal = zeros(60, 1556);
    for i = 1:60
        ZResidualsOriginal(i, :) = ZResiduals(i).d;
    end
    
    ZParametersOriginal = zeros(60, 11);
    for i = 1:60
        ZParametersOriginal(i, :) = ZParameters(i).d;
    end
    
    %% set param same as in previous implementation
    Options.radius = linspace(2,7.1,10); % radius of the cell
    Options.n_zero = linspace(1.1,1.4,10); % constant part of the refractive index
    Options.h = 0.25; % scale index in gamma range
    Options.scaleRef = true;
    Options.fixIteraionNumber = false; % Int: Fixed number of iterations
    
    tic
	[ZProcessed_test, ZResiduals_test, ZParameters_test, NumberIterations_test] = ME_EMSC(Ref, Abs, RawSpectra, wn, Options, weights);
	toc
 
    maximumArraysDifference = max( max( abs( ZProcessed_test - ZProcessedOriginal ) ) );
    
    arraysDifference = max( max( abs( ZResiduals_test - ZResidualsOriginal ) ) );
    
    if ( arraysDifference > maximumArraysDifference )
        maximumArraysDifference = arraysDifference;
    end
    
%     arraysDifference = max( max( abs( ZParameters_test ) - abs( ZParametersOriginal ) ) ); 
%     
%     if ( arraysDifference > maximumArraysDifference )
%         maximumArraysDifference = arraysDifference;
%     end
    
%     arraysDifference = max( abs( NumberIterations_test - NumberIterations ) );
%     
%     if ( arraysDifference > maximumArraysDifference )
%         maximumArraysDifference = arraysDifference;
%     end
    
    accuracy = 1e-5;
    
    if ( maximumArraysDifference > accuracy)
        assert( false )
    else
        assert( true )
    end

end