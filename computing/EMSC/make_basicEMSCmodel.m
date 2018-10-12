function [EMSCmodel] = make_basicEMSCmodel(referenceSpectra, wn, EMSCoption)
    %  Establishes the basic EMSC model without any extensions in addition to the linear and quadratic term
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
    %  referenceSpectra - Reference spectrum (row vector, or matrix containg one spectrum per row)
    %  wn               - Wavenumbers corresponding to referenceSpectra (coloumn vector)
    %  EMSCoption       - 1 (Basic EMSC, default) or 2 (MSC)
    %
    %  Output: 
    %  EMSCmodel    - Matrix containing the elements of the EMSC model as coloumn vecotrs: 
    %                 Coloumn number 1: constant baseline 
    %                 Coloumn number 2: linear
    %                 Coloumn number 3: quadratic    
    %                 Coloumn number 4: reference spectrum (mean of all spectra in referenceSpectra)


    if (nargin == 2)
        %default: all model functions
        EMSCoption = 1;
    end


    %% Calculate the basic model functions
    [~, Ny] = size(referenceSpectra);

    Start = wn(1);
    End = wn(Ny);

    C = 0.5 * ( Start + End );
    M0 = 2.0 / ( Start - End );
    M = 4.0 / ( ( Start - End ) * ( Start - End ) );

    WaveNumT = wn';
    Baseline = ones(1,Ny);
    Mean = mean(referenceSpectra,1);

    if (EMSCoption == 1) 
        Linear = M0*(Start-WaveNumT)-1;
        Quadratic = M*(WaveNumT-C).^2;
        MModel = [Baseline; Linear; Quadratic; Mean]; % EMSC: all model functions 
    elseif (EMSCoption==2)
        MModel=[Baseline;Mean]; % MSC 
    end
    
    EMSCmodel= MModel'; 
end

