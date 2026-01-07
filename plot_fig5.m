clear all, close all
load fig5data

figure('Resize', 'off', 'Units', 'centimeters', 'Position', [0 0 svgData.figureWidth svgData.figureHeight], ...
       'PaperSize', [svgData.figureWidth svgData.figureHeight], 'PaperPosition', [0 0 svgData.figureWidth svgData.figureHeight]);

population_type = 'patient';
% population_type = 'session';

set(0, 'defaultaxesfontname', 'Arial'),set(0, 'defaulttextfontname', 'Arial')
%% Compute z-values for each mere context neuron
% sn = get_unitinfo('sitenums');
% mcu = all_ps(:,2)<10^-3 & all_ps(:,1)>10^-3 & sn>0;
% l_mcu = find(mcu)
% Nm = sum(mcu);
% 
% t_min = 500;
% t_max = 1300;
% 
% resmat = [];
% blmats = [];
% qmats = [];
% im1mats = [];
% im2mats = [];
% 
% for u = 1:Nm
%     
%     s = get_unitinfo('session',all_unitinfo,l_mcu(u));
%     qs = strialinfos{s}.q;
%     disp(u)
%     bl = cellfun(@(x) sum(x>-700 & x<100), q_sptimes{l_mcu(u)})/.8;
%     M = mean(bl);
%     STD = std(bl);
%     
%     q  = cellfun(@(x) sum(x>t_min & x<t_max), q_sptimes{l_mcu(u)})/.8;
%     im1 = cellfun(@(x) sum(x>t_min & x<t_max), im1_sptimes{l_mcu(u)})/.8;
%     im2 = cellfun(@(x) sum(x>t_min & x<t_max), im2_sptimes{l_mcu(u)})/.8;
%     
%     resmat = [resmat; grpstats([(bl-M)/STD (q-M)/STD (im1-M)/STD (im2-M)/STD],qs,'mean')]; % mean z values per question for each mere context neuron
%   
%     
% end


% resmat contains z-values computed for each mere context neuron
% that were averaged per question (so each neuron contributes 5 z-values)

BLs = resmat(:,1);                      % baseline z-values
Qs = resmat(:,2);                       % question z-values
IM1s = resmat(:,3);                     % image1 z-values
IM2s = resmat(:,4);                     % image2 z-values
IMs = (resmat(:,3)+resmat(:,4))/2;      % mean image z-values
%% Compute z-values for each question for each subject/session
switch population_type
    case 'patient'
        mcu_population = mcu_subjects;
    case 'session'
        mcu_population = mcu_sessions;
end

sBLs = [];
sQs = [];
sIMs = [];

for q = 1:5 % z-values corresponding to each of the 5 questions (see resmat above)
    sBLs = [sBLs; grpstats(BLs(q:5:end),mcu_population)];
    sQs = [sQs; grpstats(Qs(q:5:end),mcu_population)];
    sIMs  = [sIMs; grpstats(IMs(q:5:end),mcu_population)];
end

[cccc cn] = get_seaborn;
%% b
% Scatter plot of mean 
% z-values for the five contextual questions during question versus picture 
% presentations, computed per MC neuron and patient


pos = svgData.rectangles.('fig5b').Position;
axes('Units', 'centimeters', 'Position', pos);


[rho,pval] = corr(sQs,sIMs)
plot(sQs,sIMs,'*','color',cccc{3})


ylim([-1.5 2]), xlim([-1.5 2])
% text(-1.4,1.75,sprintf('r = %.2f',rho),'fontsize',9)
% text(1.9,1.75,sprintf('p = %.2E',pval),'horizontalalignment','right','fontsize',9)
lsline
axis equal
grid on
switch population_type
    case 'patient'
        L = legend(sprintf('z(MC) per question per subject (N=5x%i)',size(sIMs,1)/5),'location','South','Orientation','horizontal','fontsize',7);
    case 'session'
        L = legend(sprintf('z(MC) per question per session (N=5x%i)',size(sIMs,1)/5),'location','South','Orientation','horizontal','fontsize',7);
end

XX = xlabel('z(question)_{500-1300ms}','FontSize',9);
XX.Position(2) = XX.Position(2)+.1; 
ylabel('z(pictures)_{500-1300ms}','FontSize',9)

set(gca,'fontsize',8)
%% c
% As in panel b, but baseline versus picture 
% presentations showing no correlation 

pos = svgData.rectangles.('fig5c').Position;
axes('Units', 'centimeters', 'Position', pos);

[rho_ctrl,pval_ctrl] = corr(sBLs,sIMs);
plot(sBLs,sIMs,'*','color',cccc{8})


ylim([-1.5 2]), xlim([-1.5 2])
% text(-1.4,1.75,sprintf('r = %.2f',rho_ctrl),'fontsize',9)
% text(1.9,1.75,sprintf('p = %.2E',pval_ctrl),'horizontalalignment','right','fontsize',9)
lsline
set(gca,'YAxisLocation','right')
switch population_type
    case 'patient'
        L = legend(sprintf('z(MC) per question per subject (N=5x%i)',size(sIMs,1)/5),'location','South','Orientation','horizontal','fontsize',7);
    case 'session'
        L = legend(sprintf('z(MC) per question per session (N=5x%i)',size(sIMs,1)/5),'location','South','Orientation','horizontal','fontsize',7);
end
% L.String = L.String(1:2);
axis equal
grid on

XX = xlabel('z(baseline)_{-700-100ms}','FontSize',9);
XX.Position(2) = XX.Position(2)+.1; 

ylabel('z(pictures)_{500-1300ms}','FontSize',9)
set(gca,'fontsize',8)

%% d
% Shift-corrected 
% cross-correlations of EC-MS and H-MC during picture presentations after a 
% preferred (maximum response, green) versus a non-preferred context (grey). 
pos = svgData.rectangles.('fig5d').Position;
axes('Units', 'centimeters', 'Position', pos);


YMIN = -1.5*10^-3;
YMAX = 4.5*10^-3;
p_perm = 0.05;
maxlag = 750;


tt = -maxlag:maxlag; % times of shift-corrected cross correlations


% MSEC_MCH_prefcontext: shift-corrected cross correlations between entorhinal
% mere stimulus neurons and hippocampal mere context neurons when the
% preferred context, i.e. question, of the mere context neuron 
% WAS presented (green):
h1=cplot(tt,MSEC_MCH_prefcontext,{'color',cccc{3}}); hold on,


%%
% MSEC_MCH_nonprefcontext: shift-corrected cross correlations between entorhinal
% mere stimulus neurons and hippocampal mere context neurons when the
% preferred context, i.e. question, of the mere context neuron was
% NOT presented (grey):
h2=cplot(tt,MSEC_MCH_nonprefcontext,{'color',[.3 .3 .3]});

plot_stats(tt,MSEC_MCH_pref_vs_nonpres_stats,.0008);
ylabel('cross corr. (shift-corrected)')
XX =xlabel('lag in ms');
XX.Position(2) = XX.Position(2)*1.35;
hold on 
plot([-maxlag maxlag],[0 0],'k--') 
ylim([YMIN YMAX])
text(50,YMAX-YMAX/2,sprintf('* of mere context\n   neurons in H'),'fontsize',8)

L = legend([h1.patch h2.patch], 'preferred context*','non-preferred context*','location','Northwest','fontsize',8);
L.String = L.String(1:2);
text(.91,.03,sprintf('N = %i',size(MSEC_MCH_prefcontext,1)),'horizontalalignment','center','fontsize',8,'units','normalized')

set(gca,'fontsize',8)
title('stim. (EC) -> context (H)','fontweight','normal')

%% e
% Context SVM-decoding accuracy of 
% H-MC trained with question and decoded with picture activity time locked  
% (0–100 ms) to EC-MS firing 
% HCN_context_dec_trainq_timelocked_EC <- from allDEC(:,3) averaged per
% session
HCN_context_dec_trainq_timelocked_EC = grpstats(table2array(allDEC(:,3)),pair_s);

cmap2 = lines(6);  % This provides two different colors


pos = svgData.rectangles.('fig5e').Position;
axes('Units', 'centimeters', 'Position', pos);
set(gca,'fontsize',11)


B2 = boxplot(HCN_context_dec_trainq_timelocked_EC);
ylim([.18 .2399])
hold on 
plot([0 4],[.2 .2],'k--')
set(gca,'XTick',[])
text(.55,.19,sprintf('train:        question\ndecode:   pictures\n( 0-100 ms after\nEC-MS firing )'),'FontSize',8,'Color',cmap2(5,:))
arrayfun(@(i) do_jitter(i,HCN_context_dec_trainq_timelocked_EC(:,i)),1)

% Adjust colors for each component of the boxplot
h = findobj(gca, 'Tag', 'Box');
for ii = 1
    % Select the relevant box for the ii-th column and change its color
    patch(get(h(ii), 'XData'), get(h(ii), 'YData'), cmap2(ii+4,:), 'FaceAlpha', .5);
end

ylabel('context dec. acc. (hippoc. MC)','FontSize',9)
label_ps([signrank(HCN_context_dec_trainq_timelocked_EC(:,1),.2)],max(HCN_context_dec_trainq_timelocked_EC(:,1))*1.02,1,16)

text(0.5, 0.96, sprintf('N = %i sessions',size(HCN_context_dec_trainq_timelocked_EC,1)), 'Units', 'normalized', 'HorizontalAlignment', 'center','fontsize',9);
set(gca,'fontsize',8)


%% f
% Context SVM decoding of H-MC 
% distinguishing whether preferred (red) versus non-preferred (blue) pictures of 
% the corresponding EC-MS were presented 

pos = svgData.rectangles.('fig5f').Position;
axes('Units', 'centimeters', 'Position', pos);



% allDEC(:,1:2): HCN_context_dec_ECSN_prefpicshown and HCN_context_dec_ECSN_noprefpicshown
% pair_s: corresponding sessions of each neuron pair
HCN_context_dec_ECSN_pref_nopref = grpstats(table2array(allDEC(:,1:2)),pair_s); 

% Define colormap
cmap = lines(6); 

B = boxplot(HCN_context_dec_ECSN_pref_nopref);
ylim([.15 .335])
hold on
plot([0 4],[.2 .2],'k--')
% grpstats(allDEC,pair_s)
% ylabel('context decoding (context neuron)','FontSize',12)
set(gca,'XTick',[])
text(.62,.187,sprintf('preferred picture\nof EC-MS shown'),'FontSize',8,'Color',cmap(2,:))
text(.62,.163,sprintf('non-preferred picture\nof EC-MS shown'),'FontSize',8,'Color',cmap(1,:))
arrayfun(@(i) do_jitter(i,HCN_context_dec_ECSN_pref_nopref(:,i)),1:2)

% Adjust colors for each component of the boxplot
h = findobj(gca, 'Tag', 'Box');
for ii = 1:2
    % Select the relevant box for the ii-th column and change its color
    patch(get(h(ii), 'XData'), get(h(ii), 'YData'), cmap(ii,:), 'FaceAlpha', .5);
end

pp = signrank(HCN_context_dec_ECSN_pref_nopref(:,1),HCN_context_dec_ECSN_pref_nopref(:,2),'tail','right');
Nselected = size(HCN_context_dec_ECSN_pref_nopref,1);

ss = sigstar([1,2],pp);
set(ss(2),'fontsize',16)

text(.5, 0.96, sprintf('N = %i sessions',size(HCN_context_dec_ECSN_pref_nopref,1)), 'Units', 'normalized', 'HorizontalAlignment', 'center','fontsize',9);
set(gca,'YAxisLocation','right')
ylabel('context dec. acc. (hippoc. MC)','FontSize',9)
set(gca,'fontsize',8)
%%
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%%
function plot_stats(tt,stats,yloc)
    if numel(stats.froms)>0
        plot(tt([stats.froms;stats.tos]),-[yloc yloc],'color','r','LineWidth',2)
    end
end

function do_jitter(x,ys)
    sy = size(ys);
    dataX = repmat(x,sy);
    jitterAmount = 0.05;
    jitterValuesX = 2*(rand(sy)-0.5)*jitterAmount;
    scatter(dataX+jitterValuesX,ys,10,[0 0 0],'filled','jitter','on','JitterAmount',0.06)
end