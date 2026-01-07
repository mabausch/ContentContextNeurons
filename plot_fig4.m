%%
set(0, 'defaultaxesfontname', 'Arial'),set(0, 'defaulttextfontname', 'Arial')
rng(3)
clear all, close all
load fig4data

%% plotting setup and data preparation
YMIN = -1.5*10^-3;
YMAX = 3*10^-3;

F = figure('Position',[409   323   831   655]);
subplot(2,4,1:2),subplot(2,4,3:4),subplot(2,4,5:6),SP1 = subplot(2,4,7);SP2 =subplot(2,4,8);
ax1 = get(SP1,'InnerPosition');
ax2 = get(SP2,'InnerPosition');
close(F)


% plotting colors
[cccc cn] = get_seaborn;
black = '\color[rgb]{0 0 0}';


% cross correlation times in ms
maxlag = 750;
tt = -maxlag:maxlag;


figure('Position',[409   323   831   655])

%% a left
% Mean cross-correlograms (six sessions) 
% during picture presentations between entorhinal MS and either hippocampal 
% MC neurons (blue, data_a_session) or a matched 
% number of other hippocampal neurons (non-significant, black; data_ctrl_a_session)

subplot(2,4,1:2)

h1 = cplot(tt,data_a_session,{'b-','LineWidth',2}); hold on
h2 = cplot(tt,data_ctrl_a_session,{'k-','LineWidth',2});
hold on 
plot([-maxlag maxlag],[0 0],'k--') 


% stats_a: statistics of cluster permutation test
plot_stats(tt,stats_a,.0008);
[~,iii] = max(mean(data_a_session));
text(-900,-YMIN+(YMAX-YMIN)/5,sprintf('max peak:\n%i ms',tt(iii)))
ylim([YMIN YMAX])
xlabel('lag in ms')
ylabel('cross correlation (shift-corrected)')
title('stimulus (EC) <-> context (H)','fontweight','normal')

text(.97,.03,sprintf('N = %i sessions',size(data_a_session,1)),'horizontalalignment','right','fontsize',8,'units','normalized')

%% a right
% The right panel is the same as the left panel, but 
% with the region order reversed. 
% 
subplot(2,4,3:4)
h1 = cplot(tt,data_b_session,{'b-','LineWidth',2}); hold on
h2 = cplot(tt,data_ctrl_b_session,{'k-','LineWidth',2});
hold on 
plot([-maxlag maxlag],[0 0],'k--') 

plot_stats(tt,stats_b,.0008);

ylim([YMIN YMAX])
xlabel('lag in ms')
ylabel('cross correlation (shift-corrected)')
set(gca,'YAxisLocation','right')
legend([h1.patch h2.patch],{sprintf('%smere stimulus%s vs.\n%smere context %sneurons',cn{4},black,cn{3},black),sprintf('%smere stimulus %s vs.\nother neurons',cn{4},black)})
title('stimulus (H) <-> context (EC)','fontweight','normal')

text(.03,.03,sprintf('N = %i sessions',size(data_b_session,1)),'horizontalalignment','left','fontsize',8,'units','normalized')


%% b
% Cross-correlations of all 40 
% entorhinal cortex stimulus and hippocampus context neuron pairs (MS–MC) 
% before (all_pres, grey) and after (all_posts, red) the experiment
subplot(2,4,5:6)

h1 = cplot(tt,all_pres,{'color',cccc{8}});hold on
h2 = cplot(tt,all_posts,{'color',cccc{4}});

h1.mainLine.LineWidth = 1.5;
h2.mainLine.LineWidth = 1.5;

[~,iii1] = max(mean(all_pres));
[~,iii2] = max(mean(all_posts));
text(-900,-YMIN+(YMAX-YMIN)/6.3,sprintf('max peak:\n%s%i ms\n%s%i ms',cn{8},tt(iii1),cn{4},tt(iii2)))

ttmin = -100;ttmax = -10;

plot([ttmax ttmax],[YMIN YMAX],'k--')
plot([ttmin ttmin],[YMIN YMAX],'k--')

ylabel('cross correlation (shift-corrected)')
L = legend([h1.patch h2.patch],{sprintf('%spre (%sMS%s-%sMC%s)',black,cn{4},black,cn{3},black),sprintf('%spost (%sMS%s-%sMC%s)',black,cn{4},black,cn{3},black)});
plot_stats(tt,stats_c,.0008);
ylim([[YMIN YMAX]])
xlabel('lag in ms')
title('pre versus post experiment','fontweight','normal')

text(.97,.03,sprintf('N = %i neuron pairs',size(all_pres,1)),'horizontalalignment','right','fontsize',8,'units','normalized')

%% d
% Boxplots of mean cross-correlations (−10 to 100 ms; see dashed lines) between 
% the same neuron pairs as in panel b

% Define time window mask and color indices
time_mask = (tt > ttmin) & (tt < ttmax);
color_indices = [8, 5, 1, 4];  % [pre, exp1, exp2, post]

% Compute mean cross-correlations for each experimental phase
R = [mean(all_pres(:, time_mask), 2), ...
     mean(all_exp1(:, time_mask), 2), ...
     mean(all_exp3(:, time_mask), 2), ...
     mean(all_posts(:, time_mask), 2)];

% --- Cross-correlation boxplot ---
SP1 = subplot(2, 4, 7);
YMIN = -3e-3;
YMAX = 8e-3;

B = boxplot(R);
hold on
plot([0 5], [0 0], 'k--')

style_boxplot(B, cccc, color_indices);
add_phase_labels(YMIN * 1.15, cn);

ylabel('cross correlation')
add_exp_half_legend(cn);

% Significance tests against zero
y_offset = 5.8e-3;
p_zeros = test_against_zero(R, y_offset);

% Pairwise comparisons (pre vs all others)
y_increments = (-0.2:0.1:0.5) * 1e-3;
diffmat_p = pairwise_signrank(R, 1, 1:4, y_increments);

ylim([YMIN YMAX])
set(gca, 'FontSize', 10)
text(0.14, 0.03, sprintf('N = %i', size(R, 1)), ...
    'HorizontalAlignment', 'center', 'FontSize', 8, 'Units', 'normalized')

%% Firing rates boxplot
% Firing rates before (Frs(:,1)), during (Frs(:,2:3) each half) and after (Frs(:,4))

SP2 = subplot(2, 4, 8);

B = boxplot(Frs);
style_boxplot(B, cccc, color_indices);

ylabel('firing rates (Hz)')
set(gca, 'YAxisLocation', 'right')

% Pairwise comparisons (all pairs)
y_increments = 0.1:0.15:5;
diffmat_f = pairwise_signrank(Frs, 1:4, 1:4, y_increments);

% Axis labels and limits
add_phase_labels(-2.8, cn);
ylim([-2 17.5])
add_exp_half_legend(cn);

L.String = L.String(1:2);
set(gca, 'FontSize', 10)
text(0.86, 0.03, sprintf('N = %i', size(Frs, 1)), ...
    'HorizontalAlignment', 'center', 'FontSize', 8, 'Units', 'normalized')

%% Panel labels (a, b, c)
add_panel_labels({'a', 'b', 'c'}, [0.075, 0.075, 0.525], [0.985, 0.5, 0.5]);

%% Adjust subplot positions
ax1(1) = ax1(1) + 0.03;
set(SP1, 'Position', ax1)
set(SP2, 'Position', ax2)


%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function style_boxplot(B, colors, color_indices)
% STYLE_BOXPLOT Apply custom colors and line widths to boxplot elements
    set(gca, 'XTickLabel', '')
    for i = 1:numel(color_indices)
        set(B(5, i), 'Color', colors{color_indices(i)}, 'LineWidth', 1.5)
        set(B(6, i), 'Color', colors{color_indices(i)})
    end
end


function add_phase_labels(y_pos, cn)
% ADD_PHASE_LABELS Add pre/exp/post labels below x-axis
    text(1, y_pos, sprintf('%spre', cn{8}), 'HorizontalAlignment', 'center')
    text(2.5, y_pos, 'exp.', 'HorizontalAlignment', 'center')
    text(4, y_pos, sprintf('%spost', cn{4}), 'HorizontalAlignment', 'center')
end


function add_exp_half_legend(cn)
% ADD_EXP_HALF_LEGEND Add legend for first/second experimental half
    text(0.5, 0.07, sprintf('%sfirst half', cn{5}), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'Units', 'normalized')
    text(0.5, 0.03, sprintf('%ssecond half', cn{1}), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'Units', 'normalized')
end


function p_values = test_against_zero(data, y_offset)
% TEST_AGAINST_ZERO Wilcoxon signed-rank test against zero for each column
    n_cols = size(data, 2);
    p_values = nan(1, n_cols);
    
    for j = 1:n_cols
        p_values(j) = signrank(data(:, j), 0);
        label_ps(p_values(j), y_offset - 0.1e-3, j)
        plot(j, y_offset + y_offset/15, '.w')
    end
end


function p_matrix = pairwise_signrank(data, ref_cols, compare_cols, y_increments)
% PAIRWISE_SIGNRANK Pairwise Wilcoxon signed-rank tests with sigstar annotations
    p_matrix = [];
    sigstar_handles = {};
    
    for i = ref_cols
        for j = compare_cols
            if j > i
                p = signrank(data(:, i), data(:, j));
                p_matrix = [p_matrix; i, j, p];
                if p < 0.05
                    sigstar_handles{end+1} = sigstar([i, j], p);
                end
            end
        end
    end
    
    % Adjust sigstar positions to prevent overlap
    for k = 1:numel(sigstar_handles)
        if k <= numel(y_increments)
            set(sigstar_handles{k}(1), 'YData', get(sigstar_handles{k}(1), 'YData') + y_increments(k))
            set(sigstar_handles{k}(2), 'Position', get(sigstar_handles{k}(2), 'Position') + [0, y_increments(k), 0])
        end
    end
end


function add_panel_labels(labels, x_positions, y_positions)
% ADD_PANEL_LABELS Add bold panel labels (a, b, c, etc.) to figure
    AxesH = axes('Parent', gcf, ...
        'Units', 'normalized', ...
        'Position', [0, 0, 1, 1], ...
        'Visible', 'off', ...
        'XLim', [0, 1], ...
        'YLim', [0, 1], ...
        'NextPlot', 'add');
    
    for i = 1:numel(labels)
        text(x_positions(i), y_positions(i), labels{i}, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontWeight', 'bold', ...
            'FontSize', 16, ...
            'FontName', 'SansSerif')
    end
end


function plot_stats(tt, stats, yloc)
% PLOT_STATS Plot significance bars from cluster permutation test results
    if numel(stats.froms) > 0
        plot(tt([stats.froms; stats.tos]), -[yloc yloc], 'Color', 'r', 'LineWidth', 2)
    end
end
