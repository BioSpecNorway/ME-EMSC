function mie_hulst_complex_scaledTest()

	%% Test 1: fisrt iteration input agruments check
	load 'data/mie_hulst_complex_scaled/input_dataset_1.mat';
	load 'data/mie_hulst_complex_scaled/output_dataset_1.mat';

    tic
	[ZQ_test] = mie_hulst_complex_scaled(ZRefIndexComplex.d, RefIndexWN, gammaCoef, alpha0);
	toc

    maximumArraysDifference = max( abs( ZQ_test - ZQ.d ) );
    accuracy = 1e-10;
    
    if ( maximumArraysDifference > accuracy)
        assert( false )
    else
        assert( true )
    end

end