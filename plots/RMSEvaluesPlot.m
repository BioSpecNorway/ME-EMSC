function RMSEvaluesPlot(RMSE, RMSE_limit, selectedSpectraNumbersForCorrection)
    %% Plot RMSE values
    figure; 
    stem(RMSE, 'filled'); 
    hold on 
    plot(RMSE_limit.*ones(1, length(selectedSpectraNumbersForCorrection))); 
    axis tight 
    ylim([0 max(RMSE)+0.1*max(RMSE)]);
    title('Quality test - RMSE limit') 
    set(gcf, 'Color', [1 1 1])
end

