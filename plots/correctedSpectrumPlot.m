function correctedSpectrumPlot(ProcessedQT, RefSpec, wn)
    %% Plot the corrected spectra from the last iteration step 
    figure;
    plot(wn,RefSpec,'k','Linewidth',1.5); 
    hold on 
    plot(wn, ProcessedQT', 'LineWidth', 1); 
    axis tight
    set(gcf,'Color',[1 1 1]);
    set(gca,'XDir','reverse');
    xlabel('Wavenumber [cm^-^1]','FontSize',12);
    ylabel('Absorbance','FontSize',12);
    title('Corrected spectra, reference spectrum (black)');
end

