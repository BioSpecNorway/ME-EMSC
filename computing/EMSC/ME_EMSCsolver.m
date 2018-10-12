function [Corrected, Residuals, Parameters] = ME_EMSCsolver(RawSpectra, EMSCModel)                      
    %  Solves the ME-EMSC correction. 
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
    %  RawSpectra - Raw spectra (matrix containg one spectrum per row)
    %  EMSCmodel    - Matrix containing the elements of the EMSC model as coloumn vecotrs: 
    %                 Coloumn number 1: constant baseline 
    %                 Coloumn number 2: reference spectrum 
    %                 Coloumn number 3 - end: loadings from PCA on Mie extinction curves  
    %
    %  Output: 
    %  Corrected    - Corrected spectra (matrix containg one spectrum per row)
    %  Residuals       - Residuals after correction (matrix containg one spectrum per row)
    %  Parameters     - EMSC parameters in the following order: 
    %                   constant baseline (parameter c), reference spectrum (parameter b), loadings from PCA on Mie extinction curves (parameters g1, g2 etc.)     

    Model=EMSCModel; 
    Parameters  = ( Model \ RawSpectra' )';  
    Corrected = RawSpectra - Parameters * EMSCModel' + Parameters(:,2) * EMSCModel(:,2)';
    Corrected = bsxfun(@times, Corrected, 1./Parameters(:,2)); % correct multipl. eff.
    Residuals = RawSpectra - Parameters * EMSCModel';
end