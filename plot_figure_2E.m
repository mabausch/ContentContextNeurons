%%
load all_unitinfo       % metainfo for each neuron
load rmANOVA            % data of the partial eta squares and p-values obtained by the repeated-measures ANOVA for each neuron
%%
figure('Position',[281   304   735   317])
[ celldata ] = plot_2E(all_unitinfo,all_ps,all_etas,all_etas_imperm_boot,all_etas_qperm_boot,'patient')