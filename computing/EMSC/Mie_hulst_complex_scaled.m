function [Q] = Mie_hulst_complex_scaled(RefIndexComplex, wn, gammaCoef, alphaCoef)
    % Calculates approximate extinction of electromagnetic radiation by a spehere. 
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
    %  RefIndexComplex  - Imaginary and real fluctuating part of refractive index (row vector)
    %  wn               - Wavenumbers corresponding to RefIndexComplex (coloumn vector) 
    %  gammaCoef        - Physical paramter gamma (float)
    %  alphaCoef        - Physical paramter alpha0  (float)
    %
    %  Output: 
    %  Q    - Mie extinction curve (row vector)
    
    nv = real(RefIndexComplex); 
    nvp = imag(RefIndexComplex);

    rhov = (alphaCoef * ( 1.0 + gammaCoef * nv )) .* (wn') * 100 ;  % 100 is correcting for the unit cm^-1

    divider = ( 1.0 / gammaCoef ) + nv;
    tanbeta = nvp ./ divider;
    beta0 = atan2( nvp, divider );

    cosB = cos(beta0);

    Q = 2.0 + ( 4 ./ ( rhov ) ) .* ( cosB ) .* ( - exp( -(rhov) .* (tanbeta) ) .* ( sin( (rhov) - (beta0) ) ...
       + ( 1.0 ./ ( rhov ) ) .* (cosB) .* cos( ( rhov ) - 2 .* (beta0) ) )...
       + ( 1.0 ./ ( rhov ) ) .* (cosB) .* cos( 2 .* (beta0) ) );  
end


           
       














      