function [RefSpecFitted, RawSpecFitted, wn] = adjustWavenumbers(RefSpec, wnRefSpec, RawSpec, wnRawSpec)
    %  Adjusting the wavenumbers of a reference spectrum to be compatible with the measured data set.  
    %
    %  ---------------------------------------------------------------------------------
    %  Written by: 
    %  Johanne Solheim                                     
    %                                                                                               
    %  Faculty of Science and Technology (REALTEK)                                   
    %  Norwegian Unversity of Life Sciences (www.nmbu.no)                                                                                                                             
    %                                                                                               
    %  Post address:                                                                                                                                                                         
    %  PO Box 5003, 1432 Aas, Norway                                                               
    %                                                                                                                                                               
    %  Revision: 2018  
    %  ---------------------------------------------------------------------------------
    %
    %  Input: 
    %  RefSpec      - Reference spectrum, row vector
    %  wnRefSpec    - Wavenumbers corresponding to the reference spectrum, coloumn vector 
    %  RawSpec      - Raw spectra, matrix where each row corresponds to one spectrum 
    %  wnRawSpec    - Wavenumbers corresponding to the raw spectra, coloumn vector
    %  
    %  Output: 
    %  RefSpecFitted    - Reference spectrum adjusted to the wavenumbers in wn 
    %  RawSpecFitted    - Raw spectra adjusted to the wavenumbers in wn  
    %  wn               - Wavenumbers in the range where wnRefSpec and wnRawSpec overlap, with the same spacing as wnRawSpec
   
    minWavenumber = max(min(wnRefSpec), min(wnRawSpec)); 
    maxWavenumber = min(max(wnRefSpec), max(wnRawSpec)); 

    [~,i1]=min(abs(wnRefSpec-minWavenumber));
    [~,i2]=min(abs(wnRefSpec-maxWavenumber));
    RefSpec=RefSpec(:,[i1:i2]);
    wnRefSpec = wnRefSpec(i1:i2);
    [~,j1]=min(abs(wnRawSpec-minWavenumber));
    [~,j2]=min(abs(wnRawSpec-maxWavenumber));
    RawSpec=RawSpec(:,[j1:j2]);
    wnRawSpec = wnRawSpec(j1:j2);
    
    RefSpecFitted = interp1(wnRefSpec', RefSpec, wnRawSpec');
    RawSpecFitted = RawSpec;
    wn = wnRawSpec; 
    
    if any(isnan(RefSpecFitted))
        RefSpecFitted(1)=RefSpecFitted(2);
        RefSpecFitted(end)=RefSpecFitted(end-1);
    end 
end 