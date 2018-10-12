function referenceSpectrumPlot(RefSpec, wn, options)
    % Plot of reference spectrum with weights
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
    %  RefSpec   - Initial reference spectrum, row vector 
    %  wn        - Wavenumbers corresponding to RefSpec, coloumn vector  
    %  options   - options struct with default and overwritten options 

  
    if options.Weights && strcmp(options.mode, 'Correction')
        weights = calcWeightFunction(wn, options);
    else 
        weights = ones(1,length(wn));
    end 
    
    figure;
    plot(wn, RefSpec, 'b');
    hold on 
    plot(wn, weights, 'r');
    axis tight
    ylim([-0.2 1.2]); 
    set(gcf, 'Color', [1 1 1])
    set(gca,'XDir','reverse');
    title('Reference spectrum (blue) and weight function (red)');  
    xlabel('Wavenumbers [cm^{-1}]') 
    ylabel('Absorbance') 
end