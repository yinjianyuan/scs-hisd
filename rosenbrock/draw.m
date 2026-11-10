close all;
X1={hisd_nobb(:,1),
    ahisd_nobb(:,1),
    scs_nobb(:,1),
    hisd_bb(:,1),
    ahisd_bb(:,1),
    scs_bb(:,1)}; 
Y1={hisd_nobb(:,2),
    ahisd_nobb(:,2),
    scs_nobb(:,2),
    hisd_bb(:,2),
    ahisd_bb(:,2),
    scs_bb(:,2)}; 


figure('Position', [100, 100, 700, 500], 'Color', 'w');
hold on; box on;

markers = {'o', 'd', 's', '^', 'p', 'h'};
colors = lines(6); colors = colors([1 3 5 2 4 6],:);

legend_str = {'HiSD', 'A-HiSD', 'SCS-HiSD','HiSD(BB)', 'A-HiSD(BB)', 'SCS-HiSD(BB)'};
lines_plots = [];

ii = 1:6;


for i = ii
    xx = X1{i};
    yy = Y1{i};
    
    marker_interval = 8; 
    marker_indices = 1:marker_interval:length(xx); 
    if length(xx)<12
        marker_indices = [1 5];
    end

        p = semilogy(xx, yy, ...
            'LineWidth', 1.5, ...
            'Color', colors(i,:), ...
            'Marker', markers{i}, ...
            'MarkerSize', 5, ...
            'MarkerIndices', marker_indices, ... 
            'MarkerFaceColor', 'w');
    
    lines_plots = [lines_plots, p];
end


set(gca, 'YScale', 'log');
ylim([1e-8, 10^1]); 
xlim([0, 500]); 
yticks(10.^(-10 : 2 : 0)); 
xticks(0:250:500)
set(gca, 'YMinorTick', 'off');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 1);

xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\|x^{(n)} - x^*\|_2$', 'Interpreter', 'latex', 'FontSize', 18);

legend(lines_plots, legend_str, ...
    'Interpreter', 'latex', ...
    'Location', 'northeast', ...
    'NumColumns', 2, ...
    'FontSize', 18);

hold off;

exportgraphics(gcf, ['fig' num2str(casenumber) '.pdf'], 'ContentType', 'vector');