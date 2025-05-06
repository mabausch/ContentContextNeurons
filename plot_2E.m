function [celldata] = plot_2E(all_unitinfo, all_ps, all_etas, all_etas_imperm_boot, all_etas_qperm_boot, level)
% plot_2E Generate boxplot visualization of partial eta squared values
%
% plot_2E creates a figure displaying the distribution of partial eta squared 
% values across different effect types (stimulus, context, stimulus-context),
% comparing actual data against permutation test results.
%
% Inputs:
%   all_unitinfo - Structure containing metadata for each neuron
%   all_ps - Matrix (Nunits x 3) of p-values from ANOVA tests
%   all_etas - Matrix (Nunits x 3) of partial eta squared values
%   all_etas_imperm_boot - Matrix (Nunits x 3) of eta values from image permutation
%   all_etas_qperm_boot - Matrix (Nunits x 3) of eta values from question permutation
%   level - String, level of analysis ('patient' or 'session'), default: 'patient'
%
% Outputs:
%   celldata - Cell array of data used for plotting
%
% Example:
%   celldata = plot_2E(all_unitinfo, all_ps, all_etas, all_etas_imperm_boot, 
%                    all_etas_qperm_boot, 'patient');
%
% Dependencies:
%   boxplot2.m, get_seaborn.m, get_sitenames.m, get_unitinfo.m, sigstar.m
%
% Author: Marcel Bausch
% Version: 1.0


    % Check if level is provided, otherwise use default
    if ~exist('level', 'var') || isempty(level)
        level = 'patient';
    end
    
    % Validate inputs
    if ~isstruct(all_unitinfo)
        error('plot_2E:InvalidInput', 'all_unitinfo must be a structure');
    end
    
    if ~isnumeric(all_ps) || ~isnumeric(all_etas) || ~isnumeric(all_etas_imperm_boot) || ~isnumeric(all_etas_qperm_boot)
        error('plot_2E:InvalidInput', 'all_ps, all_etas, all_etas_imperm_boot, and all_etas_qperm_boot must be numeric arrays');
    end
    
    if ~ischar(level) && ~isstring(level)
        error('plot_2E:InvalidInput', 'level must be a string');
    end
    
    % Validate dimensions
    [Nunits, numEffects] = size(all_etas);
    if size(all_ps, 1) ~= Nunits || size(all_etas_imperm_boot, 1) ~= Nunits || size(all_etas_qperm_boot, 1) ~= Nunits
        error('plot_2E:DimensionMismatch', 'Input matrices must have the same number of rows (units)');
    end
    
    % Check for required functions
    requiredFunctions = {'boxplot2', 'get_seaborn', 'get_sitenames', 'get_unitinfo', 'sigstar'};
    for i = 1:length(requiredFunctions)
        if ~exist(requiredFunctions{i}, 'file')
            error('plot_2E:MissingFunction', 'Required function %s not found in MATLAB path', requiredFunctions{i});
        end
    end
    
    % Get site information
    sitenumbers = get_unitinfo('sitenums', all_unitinfo, 1:Nunits);
    
    % Identify valid neurons (excluding session 26)
    valid_neurons = find(sitenumbers & get_unitinfo('session', all_unitinfo) ~= 26);
    fprintf('Analysis includes %d of %d total neurons\n', length(valid_neurons), Nunits);
    
    % Get recording IDs at the specified level
    rec_ids = get_unitinfo(level, all_unitinfo, valid_neurons);
    N = numel(unique(rec_ids));
    fprintf('Data grouped into %d %s-level means\n', N, level);
    
    % Initialize data storage
    ps = [];
    s_c_sc_etas = nan(3, 2, N); 
    celldata = cell(1, 6);
    k = 1;
    
    % Calculate statistics for each effect type
    for i = 1:3 % 1: stimulus effects, 2:context effects, 3: stimulus-context effects
        % Calculate means per recording ID (either subject IDs, i.e. patient IDs, or
        % session IDs
        eta = grpstats(all_etas(valid_neurons, i), rec_ids);
        
        % Select appropriate permutation test results
        switch i 
            case 1  % stimulus effects
                eta_boot = grpstats(all_etas_imperm_boot(valid_neurons, i), rec_ids);
            case {2, 3}  % context and stimulus-context effects
                eta_boot = grpstats(all_etas_qperm_boot(valid_neurons, i), rec_ids);
        end
        
        % Store results
        s_c_sc_etas(i, 1, :) = eta;
        s_c_sc_etas(i, 2, :) = eta_boot;
        
        % Perform statistical test with Bonferroni correction
        [p, ~, stats] = signrank(eta, eta_boot);
        ps = [ps 3*p]; % Bonferroni correction
        
        % Store data for output
        celldata{k} = eta;
        k = k + 1;
        celldata{k} = eta_boot;
        k = k + 1;
    end
    
    % Format p-values for display (not used in plotting, but useful for reporting)
    formatted_values = arrayfun(@format_p_value, ps, 'UniformOutput', false);
    
    % Get color maps
    [cm, cm_txt] = get_seaborn();
    
    % Create boxplot
    h = boxplot2(s_c_sc_etas, 1:3);
    
    % Configure axis
    set(gca, 'XTick', 1:3);
    
    % Set up labels with colors
    labels = {'stimulus', 'context', 'stimulus-context'};
    sel_cols = [4, 3, 1];
    for i = 1:3
        labels{i} = [cm_txt{sel_cols(i)} labels{i}];
    end
    set(gca, 'XTickLabels', labels);
    
    % Add axis labels
    ylabel(sprintf('\\eta^2 %s-means (all neurons)', level), 'FontSize', 12);
    
    % Set y-axis limits
    MAXY = max(s_c_sc_etas(:)) + 0.03;
    ylim([0 MAXY]);
    
    % Get x-values for significance markers
    x_values = arrayfun(@(x) get(x, 'XData'), h.uwhis, 'UniformOutput', 0);
    x_values = unique([x_values{:}]);
    
    % Define colors for data vs permutation tests
    cmap(1, :) = cm{5};
    cmap(2, :) = cm{8};
    
    % Apply colors to boxplot elements
    for ii = 1:2
        structfun(@(x) set(x(ii, :), 'Color', cmap(ii, :), 'MarkerEdgeColor', cmap(ii, :)), h);
    end
    
    % Add sample size text
    for ii = 1:3
        text(ii, 0.01, sprintf('(N = %i)', N), 'FontSize', 8, 'HorizontalAlignment', 'center');
    end
    
    % Configure line styles
    set([h.lwhis h.uwhis], 'LineStyle', '-');
    set(h.out, 'Marker', '.');
    
    % Add significance stars
    for i = 0:2
        x_value_pairs = x_values([1:2] + 2*i);
        a = sigstar(x_value_pairs, ps(i+1));
        hold on;
        plot([i+0.7 i+1.3], [mean(s_c_sc_etas(1, 2, :)) mean(s_c_sc_etas(1, 2, :))], '--k');
    end
    
    % Set line widths for boxplot elements
    set(h.box, 'LineWidth', 2);
    set(h.ladj, 'LineWidth', 2);
    set(h.lwhis, 'LineWidth', 2);
    set(h.med, 'LineWidth', 2);
    set(h.out, 'LineWidth', 2);
    set(h.uadj, 'LineWidth', 2);
    set(h.uwhis, 'LineWidth', 2);
    
    % Add grid
    grid on;
    
    % Add data points with jittering
    for i = 1:6
        do_jitter(x_values(i), celldata{i});
    end
    
    % Add legend
    YDOWN = 0.1;
    text(3, MAXY-0.01-YDOWN, 'data', 'HorizontalAlignment', 'center', 'Color', cmap(1, :), 'FontSize', 11);
    text(3, MAXY-0.025-YDOWN, 'other label shuffled', 'HorizontalAlignment', 'center', 'Color', cmap(2, :), 'FontSize', 11);
    text(3, MAXY-0.038-YDOWN, '(stratified)', 'HorizontalAlignment', 'center', 'Color', cmap(2, :), 'FontSize', 11);
    
    % Set font size
    set(gca, 'FontSize', 11);
    
    % Create scatter plot 1: Context vs Stimulus
    XSHIFT = 0;
    YSHIFT = -0.03;
    axes('Position', [0.44+XSHIFT 0.72-YSHIFT 0.16 0.16]);
    sel = get_unitinfo('sitenums', all_unitinfo) > 0 & get_unitinfo('session', all_unitinfo) ~= 26; % exclude session 26
    alpha_level = 0.001;
    
    plot(all_etas(sel, 2), all_etas(sel, 1), '.', 'Color', cmap(1, :));
    hold on;
    plot(all_etas_qperm_boot(sel, 2), all_etas_qperm_boot(sel, 1), '.', 'Color', cmap(2, :));
    plot([0 1], [0 1], '--k');
    
    THR = min(all_etas(all_ps(:, 2) < alpha_level & sel, 2));
    if ~isempty(THR) && ~isnan(THR)
        plot([THR THR], [0 1], 'k--', 'LineWidth', 1);
        text(THR, -0.1, '\alpha', 'HorizontalAlignment', 'left');
    end
    
    xlabel('context');
    ylabel('stimulus');
    T = title(sprintf('\\eta^2 of all %i neurons', numel(valid_neurons)), 'FontWeight', 'normal');
    
    T.Position(2) = T.Position(2) + 0.05;
    T.Position(1) = T.Position(1) + 0.9;
    ylim([0 1]);
    xlim([0 1]);
    
    % Create scatter plot 2: Stimulus-context vs Stimulus
    axes('Position', [0.70+XSHIFT 0.72-YSHIFT 0.16 0.16]);
    
    plot(all_etas(sel, 3), all_etas(sel, 1), '.', 'Color', cmap(1, :));
    hold on;
    plot([0 1], [0 1], '--k');
    plot(all_etas_qperm_boot(sel, 3), all_etas_qperm_boot(sel, 1), '.', 'Color', cmap(2, :));
    xlabel('stimulus-context');
    
    THR = min(all_etas(all_ps(:, 3) < alpha_level & sel, 3));
    if ~isempty(THR) && ~isnan(THR)
        plot([THR THR], [0 1], 'k--', 'LineWidth', 1);
        text(THR, -0.1, '\alpha', 'HorizontalAlignment', 'left');
    end
    
    ylim([0 1]);
    xlim([0 1]);
    
    % Helper function to add jittered data points
    function do_jitter(x, ys)
        sy = size(ys);
        dataX = repmat(x, sy);
        jitterAmount = 0.05;
        jitterValuesX = 2 * (rand(sy) - 0.5) * jitterAmount;
        scatter(dataX + jitterValuesX, ys, 4, [0 0 0], 'filled', 'jitter', 'on', 'JitterAmount', 0.06);
    end
    
    % Helper function to format p-values (useful for reporting)
    function formatted = format_p_value(p)
        % Format with 3 decimal places in scientific notation
        formatted = sprintf('%.3e', p);
        % Split into base and exponent parts
        parts = strsplit(formatted, 'e');
        base = str2double(parts{1});
        exp = str2double(parts{2});
        % Create final string with × symbol
        formatted = sprintf('%.3f×10%d', base, exp);
        disp(formatted);
    end
end