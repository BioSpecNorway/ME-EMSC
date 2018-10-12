function [Corrected,Residuals,Parameters] = basicEMSCsolver(RawSpectra,EMSCModel)                                                                           
    %  Solves the basic EMSC correction. 
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
    %  RawSpectra   - Raw spectra (matrix containg one spectrum per row)
    %  EMSCmodel    - Matrix containing the elements of the EMSC model as coloumn vecotrs: 
    %                 Coloumn number 1: constant baseline 
    %                 Coloumn number 2: linear
    %                 Coloumn number 3: quadratic    
    %                 Coloumn number 4: reference spectrum 
    %
    %  Output: 
    %  Corrected    - Corrected spectra (matrix containg one spectrum per row)
    %  Residuals    - Residuals after correction (matrix containg one spectrum per row)
    %  Parameters   - EMSC parameters in the following order: constant baseline, linear, quadratic, reference spectrum (parameter b)     


    Model=EMSCModel;
    Parameters  = (Model\RawSpectra')'; 
    k=1:3;
    Corrected = RawSpectra-Parameters(:,k)*EMSCModel(:,k)';
    Corrected = bsxfun(@times, Corrected, 1./Parameters(:,4)); % correct multipl. eff.
    Residuals = RawSpectra - Parameters * EMSCModel';
end