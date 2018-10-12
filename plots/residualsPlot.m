function residualsPlot(ResidualsQT, wn)
    %% Residuals at the last iteration step for every corrected spectrum
    figure; 
    set(gcf,'Color',[1 1 1]);
    plot(wn, ResidualsQT');
    set(gca,'XDir','reverse');
    axis tight
    xlabel('Wavenumber [cm^-^1]','FontSize',12);
    ylabel('Absorbance','FontSize',12);
    title('Residuals');
    hold on;
end

