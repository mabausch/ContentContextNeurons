# ContentContextNeurons

This repository contains code and data accompanying the following paper: 

### Distinct neuronal populations in the human brain combine content and context across temporal gaps

Marcel Bausch<sup>1</sup>, Johannes Niediek<sup>1</sup>, Thomas P. Reber<sup>1,2*</sup>, Sina Mackay<sup>1</sup>, Jan Boström<sup>3</sup>, Christian E. Elger<sup>1</sup>, Florian Mormann<sup>1</sup>

<sup>1</sup> Department of Epileptology, University of Bonn Medical Centre, Bonn, Germany

<sup>2</sup> Faculty of Psychology, UniDistance Suisse, Brig, Switzerland

<sup>3</sup> Department of Neurosurgery, University of Bonn Medical Centre, Bonn, Germany

If you have any questions, please contact <marcel.bausch.ukb@gmail.com>.

## Requirements

- MATLAB R2021a or later
- MATLAB Statistics and Machine Learning Toolbox

## Code Usage

To reproduce the figures in the manuscript, start MATLAB, navigate to the repository folder and add it to the MATLAB path:

```matlab
cd('/Path/to/this/repo');
addpath(pwd);

% To generate Figure 2E
plot_figure_2E

% To derive neuron counts with specific significance threshold
alpha_level = 10^-3; 
valid_neurons = find(get_unitinfo('sitenums', all_unitinfo) & get_unitinfo('session', all_unitinfo) ~= 26);
stimulus_neurons = sum(all_ps(valid_neurons, 1) < alpha_level);
context_neurons = sum(all_ps(valid_neurons, 2) < alpha_level);
stimulus_context_neurons = sum(all_ps(valid_neurons, 3) < alpha_level);

% Display results
fprintf('Neurons responding to stimulus: %d\n', stimulus_neurons);
fprintf('Neurons responding to context: %d\n', context_neurons);
fprintf('Neurons responding to stimulus-context interaction: %d\n', stimulus_context_neurons);
```

## Data Files

- **all_unitinfo.mat**: Metadata for each neuron (structure with fields 'data' and 'labels')
  - Contains information about recording sites, sessions, and patient IDs
  - Access metadata using `get_unitinfo` function, e.g., `site_info = get_unitinfo('sitenums', all_unitinfo)`

- **rmANOVA.mat**: Results from repeated-measures ANOVA for each neuron
  - **all_etas** (Nunits × 3): Partial eta squared values for stimulus, context, and stimulus-context effects
  - **all_etas_imperm_boot** (Nunits × 3): Eta values from permuting image labels PER question (control)
  - **all_etas_qperm_boot** (Nunits × 3): Eta values from permuting question labels PER image (control)
  - **all_ps** (Nunits × 3): P-values for stimulus, context, and stimulus-context effects

## Main Functions

- **plot_2E.m**: Generates boxplots comparing partial eta squared values across effect types
- **boxplot2.m**: Enhanced boxplot with better formatting than MATLAB's built-in function
- **get_seaborn.m**: Creates color palettes similar to Python's Seaborn
- **get_sitenames.m**: Returns site name identifiers
- **get_unitinfo.m**: Retrieves metadata about specific units
- **sigstar.m**: Adds significance stars to plots

## Analysis Details

### Neuron Classification

Neurons are classified as responsive to different effects based on significance thresholds:

```matlab
% Example classification with p < 0.001
alpha_level = 10^-3;
valid_neurons = find(get_unitinfo('sitenums', all_unitinfo) & get_unitinfo('session', all_unitinfo) ~= 26);

% Count neurons responding to each effect
stimulus_neurons = sum(all_ps(valid_neurons, 1) < alpha_level);
context_neurons = sum(all_ps(valid_neurons, 2) < alpha_level);
stimulus_context_neurons = sum(all_ps(valid_neurons, 3) < alpha_level);

% Count neurons responding to combinations of effects
stimulus_only = sum(all_ps(valid_neurons, 1) < alpha_level & all_ps(valid_neurons, 2) >= alpha_level & all_ps(valid_neurons, 3) >= alpha_level);
context_only = sum(all_ps(valid_neurons, 1) >= alpha_level & all_ps(valid_neurons, 2) < alpha_level & all_ps(valid_neurons, 3) >= alpha_level);
interaction_only = sum(all_ps(valid_neurons, 1) >= alpha_level & all_ps(valid_neurons, 2) >= alpha_level & all_ps(valid_neurons, 3) < alpha_level);
```

### Statistical Analysis

The analysis uses repeated-measures ANOVA to calculate partial eta squared values and p-values for:
1. Main effect of stimulus (image identity)
2. Main effect of context (question type)
3. Interaction effect between stimulus and context

Statistical significance in Figure 2E is assessed using Wilcoxon signed-rank tests with Bonferroni correction, comparing actual data against permutation test controls.