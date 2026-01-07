%% c
% A heat map of context decoding 
% across time. The green diagonal line denotes identical training or testing times 
% (see panel d).
clear all
load fig3Cdata.mat

disp(t(steps));

figure
imagesc(mean_context_decmat*100,[20 30]) 

c = colorbar();
ccc = parula(5);

hold on, plot(steps([1 1]),[0 67],'w--','LineWidth',2)
plot(steps([3 3]),[0 67],'w--','LineWidth',2)
plot(steps([5 5]),[0 67],'w--','LineWidth',2)
plot([0 67],steps([1 1]),'w--','LineWidth',2)
plot([0 67],steps([3 3]),'w--','LineWidth',2)
plot([0 67],steps([5 5]),'w--','LineWidth',2)


plot(steps([4 4]),[0 67],'y--','LineWidth',1)
plot(steps([6 6]),[0 67],'y--','LineWidth',1)
plot([0 67],steps([4 4]),'y--','LineWidth',1)
plot([0 67],steps([6 6]),'y--','LineWidth',1)


plot([0 69],[0 69],'-','LineWidth',2,'Color',ccc(3,:))


text(steps(1)+2,4,'question','color',[1 1 1])
text(steps(3)+2,4,'pic 1','color',[1 1 1])
text(steps(5)+2,4,'pic 2','color',[1 1 1])

t1 = text(4,steps(1)+2,'question','color',[1 1 1]);
t2 = text(4,steps(3)+2,'pic 1','color',[1 1 1]);
t3 = text(4,steps(5)+2,'pic 2','color',[1 1 1]);
set(t1,'Rotation',90);set(t2,'Rotation',90);set(t3,'Rotation',90);

set(gca,'XTick',steps); set(gca,'XTicklabel',{'q','1000','p1','1000','p2','1000',''},'fontsize',12)
set(gca,'YTick',[]);

ylabel(c, 'context decoding accuracy (%)','fontsize',14)
ylabel('t(trained)','fontsize',14)
xlabel('t(decoded) in ms','fontsize',14)
set(gca,'YDir','normal')


rectangle('Position',[11 12 4 4],'LineWidth',2,'EdgeColor',[.4 .4 .4])
rectangle('Position',[28 12 4 4],'LineWidth',2,'EdgeColor',[133 63 63]/255)
rectangle('Position',[44 12 4 4],'LineWidth',2,'EdgeColor',[64 67 130]/255)


title(sprintf('session-wise decoding (\\color[rgb]{0.3333 0.6588 0.4078}context neurons\\color[rgb]{0 0 0})'),'fontweight','normal','fontsize',11)

%% d
% Patient-wise context (top; green, patwise_diag_contextdec) or stimulus decoding (bottom; picture 1 or picture 2 in grey or lavender, 
% respectively, patwise_diag_im1dec and patwise_diag_im2dec) and label-shuffled controls (red, patwise_diag_labelshuffle_ctrl). 
clear all
load fig3Ddata
load fig3Ddata_extended
figure('Position',[680   558   513   420])

axes('Position',[.06 .4 .8 .5])

ccc = parula(5);
c1 = ccc(3,:);
c2 = [214, 95, 95]/255;

N = numel(t);

H = cplot(1:N,patwise_diag_contextdec,{'color',c1,});hold on, 
H2 = cplot(1:N,patwise_diag_labelshuffle_ctrl,{'color',c2});

set(H.patch,'facealpha',.4,'facecolor',c1)
set(H2.patch,'facealpha',.4,'facecolor',c2)

M = mean(patwise_diag_contextdec);

ymax = max(M);
YMIN = .16;
YMAX = ymax+.06;

plot(steps([1 1]),[YMIN YMAX],'--','LineWidth',1,'Color',[.4 .4 .4])
plot(steps([3 3]),[YMIN YMAX],'--','LineWidth',1,'Color',[.4 .4 .4])
plot(steps([5 5]),[YMIN YMAX],'--','LineWidth',1,'Color',[.4 .4 .4])

chance_H = plot([0 N],[.2 .2],'b--','LineWidth',2,'color',[.1 .1 .1]);


% fig3D_q_vs_shuflle_stats: permttest stats of context decoding vs shuffle decoding (0.25)

for u = 1:numel(fig3D_q_vs_shuflle_stats.froms)
    H_significant = plot([(fig3D_q_vs_shuflle_stats.froms(u)) (fig3D_q_vs_shuflle_stats.tos(u))],[.18 .18],'LineWidth',2,'Color',c1);
end
xlim([0 N])
xlim([0 N-2])

set(gca,'YAxisLocation','right')

YYY = ylabel('decoding accuracy','fontsize',12);
YYY.Units = 'norm';
YYY.Position(2) = .2;

[hleg1, hobj1] = legend([H.patch,H2.patch,chance_H],{'context','label shuffle','chance level'},'location','Northwest','fontsize',11);


set(gca,'XTick',steps); set(gca,'XTicklabel',{'q','1000','p1','1000','p2','1000',''})%t(steps));
set(gca,'fontsize',10)
grid on
ylim([YMIN YMAX])
cccc = get_seaborn;

text(N-12,.325,sprintf('\\color[rgb]{0.3333 0.6588 0.4078}context neurons\n\\color[rgb]{0 0 0}(%i subjects)',size(patwise_diag_contextdec,1)),'horizontalalignment','center','fontsize',10)
title('subject-wise decoding','fontweight','normal','fontsize',11)


axes('Position',[.06 .1 .8 .3])
H3 = cplot(1:N,patwise_diag_im1dec,{'color',cccc{8}});hold on
H2 = cplot(1:N,patwise_diag_im2dec,{'color',cccc{5}});


M = mean(patwise_diag_im1dec);
ymax = max(M);
YMIN = .15; 
YMAX = ymax+.1;
plot(steps([1 1]),[YMIN YMAX],'--','LineWidth',1,'Color',[.4 .4 .4])
plot(steps([3 3]),[YMIN YMAX],'--','LineWidth',1,'Color',[.4 .4 .4])
plot(steps([5 5]),[YMIN YMAX],'--','LineWidth',1,'Color',[.4 .4 .4])

chance_H = plot([0 N],[.25 .25],'b--','LineWidth',2,'color',[.1 .1 .1]);

% fig3Dim1stats: permttest stats of im1 decoding vs chance decoding (0.25)
% fig3Dim2stats: permttest stats of im2 decoding vs chance decoding  (0.25)
for u = 1:numel(fig3Dim1stats.froms)
    H_significant = plot([(fig3Dim1stats.froms(u)) (fig3Dim1stats.tos(u))],[YMIN YMIN]+.05,'LineWidth',2,'Color',cccc{8});
end
for u = 1:numel(fig3Dim2stats.froms)
    H_significant = plot([(fig3Dim2stats.froms(u)) (fig3Dim2stats.tos(u))],[YMIN YMIN]+.03,'LineWidth',2,'Color',cccc{5});
end
xlim([0 N])


set(gca,'YAxisLocation','right')

[hleg1, hobj1] = legend([H3.patch,H2.patch],{'pic 1','pic 2'},'location','Northwest','fontsize',11);


set(gca,'XTick',steps); set(gca,'XTicklabel',{'q','1000','p1','1000','p2','1000',''})%t(steps));
set(gca,'fontsize',10)

grid on
ylim([YMIN YMAX])
xlim([0 N-2])


YYY.FontSize = 13;

text(N-12,.55,sprintf('\\color[rgb]{0.7686 0.3059 0.3216}stimulus neurons\n\\color[rgb]{0 0 0}(%i subjects)',size(patwise_diag_im1dec,1)),'horizontalalignment','center','fontsize',10)
