function [ExtinctionCurves]=make_MieExtinction_curves(absorbanceSpectrum, wn, alpha0, gamma)
    %  Calculate Mie extinction curves from a given imaginary part of the complex refractive index 
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
    %  absorbanceSpectrum   - Absorbance spectrum, input for the imaginary part of the refractive index (row vector)
    %  wn                   - Wavenumbers corresponding to absorbanceSpectrum (coloumn vector) 
    %  alpha0               - Physical paramter alpha0 (range of values, row vector)
    %  gamma                - Physical paramter gamma (range of values, row vector)
    %
    %  Output: 
    %  ExtinctionCurves     - Matrix containing the Mie extinction curves (matrix, each row corresponds to one curve)
    
    % Imaginary part of the refractive index 
    RefIndexABS = absorbanceSpectrum';
    
    % Number of points for extending the absorbance spectrum 
    xts = 200; 
    

    %% Calculate the refractive index by Hilbert transform
    % Extend the absorbance spectrum in both directions 
    dx = wn(2) - wn(1); 
    wn = [(dx.*(linspace(1,xts,xts))+(wn(1)-dx*(xts+1)))'; wn; (dx.*(linspace(1,xts,xts))+wn(end))']; 
    RefIndexABS = [repmat(RefIndexABS(1), [xts 1]); RefIndexABS; repmat(RefIndexABS(end), [xts 1])]; 
    
    % Calculate real fluctuating part 
    RefIndexN = kkre_hilbert(RefIndexABS ./ ( wn * 100 ) );  

    RefIndexN = RefIndexN((xts+1):end-xts); 
    RefIndexABS = RefIndexABS((xts+1):end-xts); 
    wn = wn((xts+1):end-xts); 

    if  abs( min( RefIndexN ) ) > 1
        RefIndexIN = RefIndexN / abs( min( RefIndexN ) );  % Real fluctuating part normalized
        RefIndexABS = RefIndexABS ./ ( wn * 100 ) / abs( min( RefIndexN ) ); % Imaginary part normalized 
    else
        RefIndexIN = RefIndexN;  
        RefIndexABS = RefIndexABS ./ ( wn * 100 ); 
    end

    b = complex(0,1);

    refIndexComplex = RefIndexIN' + b * RefIndexABS';

    ifunc = 0;

    ExtinctionCurves = zeros( size(gamma, 2) * size(alpha0, 2), size(wn, 1) );

    for i=1:(size(gamma,2))
        for j=1:(size(alpha0,2))
            ifunc = ifunc + 1;
            ExtinctionCurves(ifunc,:) = Mie_hulst_complex_scaled( refIndexComplex, wn, gamma(j), alpha0(i) );  
        end
    end

end


