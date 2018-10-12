function [MieEMSCmodel, PCnumber] = make_ME_EMSCmodel(referenceSpectrum, pureAbsorbanceSpectrum, wn, alpha0, gamma, options, PCnumber)
    % Establish the ME-EMSC model. 
    %
    %  ---------------------------------------------------------------------------------
    %  Written by: 
    %  Johanne Solheim, Evgeniy Gunko, Tatiana Konevskikh, Achim Kohler                                             
    %                                                                                               
    %  Faculty of Science and Technology (REALTEK)                                   
    %  Norwegian Unversity of Life Sciences (www.nmbu.no)                                                                                                                             
    %                                                                                               
    %  Post address:                                                                                                                                                                         
    %  PO Box 5003, 1432 Aas, Norway                                                               
    %                                                                                               
    %  ---------------------------------------------------------------------------------
    %
    %  Input: 
    %  referenceSpectrum        - Reference spectrum for ME-EMSC model (row vector) 
    %  pureAbsorbanceSpectrum   - Input spectrum for calculation of n' (row vector)  
    %  wn                       - Wavenumbers corresponding to referenceSpectrum and pureAbsorbanceSpectrum (coloumn vector) 
    %  alpha0                   - Physical paramter alpha0 (range of values, row vector)
    %  gamma                    - Physical paramter gamma (range of values, row vector)
    %  options                  - options for the correction, specified in the ME_EMSC function (struct)
    %  PCnumber                 - Number of principal components in the Mie meta-model (int) 
    %
    %  Output: 
    %  MieEMSCmodel     - Matrix containing the elements of the ME-EMSC as coloumn vecotrs 
    %                     Coulmn number 1: baseline 
    %                     Column number 2: reference spectrum 
    %                     Column number 3 - end: loadings from PCA on Mie extinction curves 
    %  PCnumber         - Number of principal components in the Mie meta-model, if not specified as input (int)

	%% Create meta model (100 extinction curves)
	[MieExtinctionCurves] = make_MieExtinction_curves(pureAbsorbanceSpectrum, wn, alpha0, gamma); 

	%% Remove the mean-centered reference spectrum from the model spectra
    m = referenceSpectrum * referenceSpectrum';
	norm = sqrt( m ) ;
	rnorm = referenceSpectrum / norm;
	s = MieExtinctionCurves * rnorm' ;
	MieExtinctionCurves = MieExtinctionCurves - s * rnorm;

	%% decompose the set of Mie functions
	[lds, ~,latent] = pca(MieExtinctionCurves, 'Centered', false); 
    
    %% Construct the Mie model with the new ref spec
	[EMSCfunctions] = make_basicEMSCmodel(referenceSpectrum, wn, 2);

    if nargin < 7     
        t = zeros(1, size(latent,1)-1);
        
        for i=1:size(latent,1)-1
            s(i) = sum(latent(1:i,1));   
            t(i) = ( s(i) / sum(latent(:,1)) ) * 100;
        end

        if options.PCnumber
            PCnumber = options.PCnumber; 
        elseif options.ExplainedVariance 
            [~,id2]=find( t > options.ExplainedVariance );
            PCnumber=min(id2);
        else
            msg = 'Set options.PCnumber or options.ExplainedVariance'; 
            error(msg) 
        end 
    end
    
    [nRow, nCol] = size(EMSCfunctions);
    PCzeros = zeros( nRow, PCnumber);
    MieEMSCmodel = [EMSCfunctions, PCzeros];
	for i=1:PCnumber
	    MieEMSCmodel(:, nCol + i )= lds(:,i)';  
    end 
end