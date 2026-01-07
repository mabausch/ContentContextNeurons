function plot3A()
% PLOT3A Generate Figure 3A: Population SVM-decoding accuracies of context
%
% This function creates boxplots showing context decoding accuracies during
% picture presentations, comparing different contextual questions using
% all neurons versus only context neurons.
%
%
% Required data file: fig3Adata.mat (containing population_cell, patids)
% Required functions: boxplot2, get_seaborn, sigstar, label_ps, grpstats


%% ========================================================================
%  CONFIGURATION
%  ========================================================================

% Set default fonts for all figures
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

% Analysis parameters
CONFIG = struct();
CONFIG.populations       = {'all', 'context'};  % Neuron populations to analyze
CONFIG.analysis_window   = '100-1000';          % Time window in ms (not used in plot, for reference)
CONFIG.aggregate_by_patient = true;             % If true, average sessions per patient
CONFIG.exclude_session_26   = true;             % Exclude session 26 (behavioral exclusion)

% Dataset dimensions
CONFIG.n_patients  = 17;   % Total patients (one excluded due to behavior)
CONFIG.n_sessions  = 50;   % Total recording sessions
CONFIG.n_questions = 5;    % Number of contextual questions

% Statistical thresholds
CONFIG.chance_level        = 0.2;   % Chance level for 5-way classification (1/5)
CONFIG.significance_alpha  = 0.05;  % Alpha for between-group comparisons

% Visual parameters
VISUAL = struct();
VISUAL.figure_position = [520, 410, 774, 388];  % [left, bottom, width, height]
VISUAL.boxplot_linewidth = 2;
VISUAL.axis_fontsize     = 12;
VISUAL.ylabel_fontsize   = 11;
VISUAL.legend_fontsize   = 9;
VISUAL.annotation_fontsize = 9;

% Y-axis limits and offsets (differ by aggregation mode)
if CONFIG.aggregate_by_patient
    VISUAL.sigstar_y_offset    = 0.17;  % Vertical offset for between-group significance bars
    VISUAL.chance_test_y_pos   = 0.55;  % Y position for chance-level significance stars
else
    VISUAL.sigstar_y_offset    = 0.07;
    VISUAL.chance_test_y_pos   = 0.65;
end

% Question labels for display
QUESTION_LABELS = {'Big?', 'Last seen?', 'Older? / More expensive?', 'Like better?', 'Brighter?'};
QUESTION_LABELS_WITH_ALL = [QUESTION_LABELS, {'all'}];

%% ========================================================================
%  LOAD DATA
%  ========================================================================

% Load session-wise context decoding accuracies for each of the 5 questions
% population_cell: 1x2 cell array, each containing 50x5 matrix of decoding accuracies
% patids: session-to-patient mapping vector
load('fig3Adata.mat', 'population_cell', 'patids');

%% ========================================================================
%  DATA PREPROCESSING
%  ========================================================================

n_populations = numel(CONFIG.populations);

% Exclude session 26 if specified (behavioral exclusion criterion)
if CONFIG.exclude_session_26
    for pop_idx = 1:n_populations
        population_cell{pop_idx}(26, :) = NaN;
    end
end

% Always compute patient-level averages (needed for significance tests)
all_neurons_by_patient     = grpstats(population_cell{1}, patids);
context_neurons_by_patient = grpstats(population_cell{2}, patids);

% Aggregate data either by patient or keep session-level
if CONFIG.aggregate_by_patient
    % Initialize matrix: [population x (questions + mean) x patients]
    decoding_matrix = nan(n_populations, CONFIG.n_questions + 1, CONFIG.n_patients);
    
    % Fill in question-specific accuracies (patient-averaged)
    decoding_matrix(1, 1:CONFIG.n_questions, :) = all_neurons_by_patient';
    decoding_matrix(2, 1:CONFIG.n_questions, :) = context_neurons_by_patient';
    
    % Add mean across all questions as 6th entry
    decoding_matrix(1, 6, :) = mean(all_neurons_by_patient, 2)';
    decoding_matrix(2, 6, :) = mean(context_neurons_by_patient, 2)';
else
    % Keep session-level data
    decoding_matrix = nan(n_populations, CONFIG.n_questions + 1, CONFIG.n_sessions);
    
    decoding_matrix(1, 1:CONFIG.n_questions, :) = population_cell{1}';
    decoding_matrix(2, 1:CONFIG.n_questions, :) = population_cell{2}';
    
    decoding_matrix(1, 6, :) = mean(population_cell{1}, 2)';
    decoding_matrix(2, 6, :) = mean(population_cell{2}, 2)';
end

%% ========================================================================
%  DEFINE COLORS
%  ========================================================================

% Question-specific colors (Seaborn-inspired palette)
question_colors = {
    [0.7686, 0.3059, 0.3216], ...  % Big? (reddish)
    [0.3333, 0.6588, 0.4078], ...  % Last seen? (greenish)
    [0.2980, 0.4471, 0.6902], ...  % Older/More expensive? (bluish)
    [0.7714, 0.6753, 0.5106]       % Like better? (brownish)
};

% Get additional colors from Seaborn palette
[seaborn_colors, seaborn_color_strings] = get_seaborn();
question_colors{5} = seaborn_colors{6};  % Brighter? color

% Colors for annotations
COLORS = struct();
COLORS.mean_box     = [0.6, 0.6, 0.6];  % Gray for "all" boxplots
COLORS.mean_label   = [0.4, 0.4, 0.4];  % Darker gray for "all" text labels

%% ========================================================================
%  CREATE FIGURE
%  ========================================================================

close all;
fig = figure('Position', VISUAL.figure_position);

% Create boxplots
boxplot_handles = boxplot2(decoding_matrix, 1:n_populations);

%% ========================================================================
%  CONFIGURE AXES AND LABELS
%  ========================================================================

% Set x-axis tick labels with color formatting for "context"
population_labels = CONFIG.populations;
x_tick_labels = {
    sprintf('%s neurons', population_labels{1}), ...
    sprintf('%s%s neurons', seaborn_color_strings{3}, population_labels{2})
};

set(gca, 'XTick', 1:n_populations);
set(gca, 'XTickLabel', x_tick_labels);

ylabel('context decoding accuracy', 'FontSize', VISUAL.ylabel_fontsize);
ylim([-0.1, 1.2]);
xlim([0.2, n_populations + 0.8]);

% Draw chance level line
hold on;
plot([0, 4], [CONFIG.chance_level, CONFIG.chance_level], 'k--');
box on;

%% ========================================================================
%  APPLY COLORS TO BOXPLOTS
%  ========================================================================

% Apply question-specific colors to each boxplot group
for q_idx = 1:CONFIG.n_questions
    structfun(@(x) set(x(q_idx, :), ...
        'Color', question_colors{q_idx}, ...
        'MarkerEdgeColor', question_colors{q_idx}), boxplot_handles);
end

% Apply gray color to the "all questions" boxplots (6th group)
structfun(@(x) set(x(6, :), ...
    'Color', COLORS.mean_box, ...
    'MarkerEdgeColor', COLORS.mean_box), boxplot_handles);

%% ========================================================================
%  EXTRACT X-POSITIONS FOR ANNOTATIONS
%  ========================================================================

% Get x-coordinates of each boxplot from whisker positions
x_coords_cell = arrayfun(@(x) get(x, 'XData'), boxplot_handles.uwhis, 'UniformOutput', false);
x_positions = unique([x_coords_cell{:}]);

% Add "all" labels below the mean boxplots
text(x_positions(6), 0.15, 'all', ...
    'HorizontalAlignment', 'center', 'Color', COLORS.mean_label, 'FontWeight', 'bold');
text(x_positions(12), 0.15, 'all', ...
    'HorizontalAlignment', 'center', 'Color', COLORS.mean_label, 'FontWeight', 'bold');

%% ========================================================================
%  FORMAT BOXPLOT ELEMENTS
%  ========================================================================

% Set whisker line style to solid
set([boxplot_handles.lwhis, boxplot_handles.uwhis], 'LineStyle', '-');

% Set outlier marker style
set(boxplot_handles.out, 'Marker', '.');

% Apply consistent line width to all boxplot elements
boxplot_elements = {'box', 'ladj', 'lwhis', 'med', 'out', 'uadj', 'uwhis'};
for elem = boxplot_elements
    set(boxplot_handles.(elem{1}), 'LineWidth', VISUAL.boxplot_linewidth);
end

%% ========================================================================
%  BETWEEN-GROUP SIGNIFICANCE TESTS (Mann-Whitney U)
%  ========================================================================

% Collect handles for repositioning significance markers
sigstar_line_handles = [];
sigstar_text_handles = [];

% Compare decoding accuracies between questions within each population
for pop_idx = 1:n_populations
    for q1 = 1:CONFIG.n_questions
        % Get data for first question
        if CONFIG.aggregate_by_patient
            data1 = grpstats(population_cell{pop_idx}(:, q1), patids);
        else
            data1 = population_cell{pop_idx}(:, q1);
        end
        data1 = data1(~isnan(data1));
        
        % Compare with each other question (in specific order for layout)
        for q2 = [4, 3, 2, 1, 5]
            if q1 > q2
                % Get data for second question
                if CONFIG.aggregate_by_patient
                    data2 = grpstats(population_cell{pop_idx}(:, q2), patids);
                else
                    data2 = population_cell{pop_idx}(:, q2);
                end
                data2 = data2(~isnan(data2));
                
                % Perform Mann-Whitney U test (ranksum)
                p_value = ranksum(data1, data2);
                
                % Add significance markers if significant
                if p_value < CONFIG.significance_alpha
                    if pop_idx == 1
                        sigstar_handles = sigstar([x_positions(q1), x_positions(q2)], p_value, 0);
                    else
                        sigstar_handles = sigstar([x_positions(q1 + 6), x_positions(q2 + 6)], p_value, 0);
                    end
                    sigstar_line_handles = [sigstar_line_handles, sigstar_handles(1)];
                    sigstar_text_handles = [sigstar_text_handles, sigstar_handles(2)];
                end
            end
        end
    end
end

% Adjust vertical position of significance markers
arrayfun(@(h) set(h, 'YData', get(h, 'YData') + VISUAL.sigstar_y_offset), sigstar_line_handles);
arrayfun(@(h) set(h, 'Position', get(h, 'Position') + [0, VISUAL.sigstar_y_offset, 0]), sigstar_text_handles);

%% ========================================================================
%  CHANCE-LEVEL SIGNIFICANCE TESTS (Wilcoxon signed-rank)
%  ========================================================================

all_signrank_pvalues = [];
all_decoding_means = [];
x_pos_idx = 1;

for pop_idx = 1:n_populations
    % Test each question against chance level
    for q_idx = 1:CONFIG.n_questions
        % Get data
        if CONFIG.aggregate_by_patient
            data = grpstats(population_cell{pop_idx}(:, q_idx), patids);
        else
            data = population_cell{pop_idx}(:, q_idx);
        end
        data = data(~isnan(data));
        
        % Wilcoxon signed-rank test against chance level
        p_value = signrank(data, CONFIG.chance_level);
        
        % Store results
        all_decoding_means = [all_decoding_means, nanmean(data)];
        all_signrank_pvalues = [all_signrank_pvalues, p_value];
        
        % Add significance marker and print to console
        label_ps(p_value, VISUAL.chance_test_y_pos, x_positions(x_pos_idx));
        fprintf('%s: %.3E\n', QUESTION_LABELS_WITH_ALL{q_idx}, p_value);
        x_pos_idx = x_pos_idx + 1;
    end
    
    % Test mean across all questions against chance
    if pop_idx == 1
        p_all = signrank(mean(all_neurons_by_patient, 2), CONFIG.chance_level);
    else
        p_all = signrank(mean(context_neurons_by_patient, 2), CONFIG.chance_level);
    end
    
    label_ps(p_all, VISUAL.chance_test_y_pos, x_positions(x_pos_idx));
    fprintf('%s: %.3E\n', QUESTION_LABELS_WITH_ALL{6}, p_all);
    
    % Add sample size annotation below x-axis
    n_samples = sum(~isnan(data));
    if CONFIG.aggregate_by_patient
        sample_label = sprintf('%i subjects', n_samples);
    else
        sample_label = sprintf('%i sessions', n_samples);
    end
    text(pop_idx, -0.055, sample_label, ...
        'HorizontalAlignment', 'center', 'FontSize', VISUAL.annotation_fontsize);
    
    x_pos_idx = x_pos_idx + 1;
end

%% ========================================================================
%  FINAL FORMATTING
%  ========================================================================

% Create legend for question colors
legend_handle = legend(boxplot_handles.box(1:5), QUESTION_LABELS, ...
    'Location', 'north', 'FontSize', VISUAL.legend_fontsize, 'Orientation', 'horizontal');

% Fix legend labels for MATLAB R2021 compatibility
matlab_version = version;
if contains(matlab_version, 'R2021')
    legend_handle.String = QUESTION_LABELS;
end

% Fine-tune legend position
legend_handle.Position(1) = legend_handle.Position(1) + 0.007;

% Set axis font size
set(gca, 'FontSize', VISUAL.axis_fontsize);
box on;

% Add chance level annotation
text(2.57, 0.17, 'chance');

end