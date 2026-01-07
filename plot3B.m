function plot3B()
    %% b
    % Pooled context (green) or 
    % stimulus (red) decoding accuracies from context (left) and stimulus (right) 
    % neurons, for cross-validation or generalization across pictures (red), contexts 
    % (green) and serial picture positions (blue). 
    
    load fig3Bdata.mat
    YP=0.23;
    figure('Position',[680   664   227+154   314])
    
    axes('units','pixel','Position',[30.5100+22   35.5400  153.5917-12  255.9100])
    Nrep = 30;
    [sb,ss] = get_seaborn;
    
    % 30 subsamples, pooled decoding of context either via cross-val,
    % across pictures or across positions
    decmat1 = [mean(contextdec,2) mean(contextdec_imgen,2)  (mean(contextdec_posgen_im1,2)+mean(contextdec_posgen_im2,2))/2]; 
    
    % 30 subsamples, pooled decoding of stimulus either via cross-val,
    % across contexts or across positions
    decmat2 = [mean(stimulusdec,2) mean(stimulusdec_contextgen,2)  (mean(stimulusdec_posgen_im1,2)+mean(stimulusdec_posgen_im2,2))/2];
    
    B = boxplot(decmat1);ylim([-.1 1.2]);
    set(gca,'Fontsize',10)
    %%
    set(gca,'YTicklabel',[])
    
    hold on
    plot([1.5 1.5],[.25 1.1],'--','color',[.4 .4 .4])
    plot([2.5 2.5],[.25 .28],'--','color',[.6 .6 .6])
    % plot([1.5 1.5],[.4 1.1],'--','color',[.4 .4 .4])
    plot([2.5 2.5],[.4 1.1],'--','color',[.6 .6 .6])
    plot([0 4],[.2 .2],'k--','Linewidth',1)
    set(gca,'XTick',[])
    Fs = 9.5;
    ps = arrayfun(@(i) signrank(decmat1(:,i),.2),1:3);
    label_ps(ps,1.03)
    
    text(2.5,.17+YP,sprintf('generalization'),'horizontalalignment','center','fontsize',Fs-1,'verticalalignment','top')
    text(2.5,.12+YP,'across','horizontalalignment','center','fontsize',Fs-1,'verticalalignment','top')
    text(2,.07+YP,'pictures','horizontalalignment','center','fontsize',Fs-1,'color',sb{4},'verticalalignment','top')
    text(3,.07+YP,'positions','horizontalalignment','center','fontsize',Fs-1,'color',sb{1},'verticalalignment','top')
    
    

    
    arrayfun(@(x) set(x,'Color',sb{3}),B(:,1))
    arrayfun(@(x) set(x,'Color',sb{4}),B(:,2))
    arrayfun(@(x) set(x,'Color',sb{1}),B(:,3))
    
    
    text(.58,.16,'chance','fontsize',8)
    xlim([.5 3.5])
    text(2,1.15,'\color[rgb]{0.3333 0.6588 0.4078}context decoding','horizontalalignment','center','fontsize',9)

    text(2,.0,[sprintf('30 subsamples of\n%i%s context neurons\\color[rgb]{0 0 0} ',Nselunits,ss{3}) '(75%)'],'horizontalalignment','center','fontsize',8)
    
    
    AX2_POS = [30.5100+158+6   35.5400  153.5917-12  255.9100];
    axes('units','pixel','Position',AX2_POS)
    B = boxplot(decmat2);ylim([-.1 1.2]);
    set(gca,'Fontsize',9)
    ylabel('decoding accuracy','fontsize',11) %
    hold on
    plot([1.5 1.5],[.25+0.05 1.1],'--','color',[.4 .4 .4])
    plot([2.5 2.5],[.25+0.05 .28+0.05],'--','color',[.6 .6 .6])
    plot([2.5 2.5],[.4+0.06 1.1],'--','color',[.6 .6 .6])
    plot([0 4],[.25 .25],'k--','Linewidth',1)
    set(gca,'XTick',[])
    
    Fs = 9.5;
    text(2.5,.22+YP,sprintf('generalization'),'horizontalalignment','center','fontsize',Fs-1,'verticalalignment','top')
    text(2.5,.17+YP,'across','horizontalalignment','center','fontsize',Fs-1,'verticalalignment','top')
    text(2,.12+YP,'contexts','horizontalalignment','center','fontsize',Fs-1,'color',sb{3},'verticalalignment','top')
    text(3,.12+YP,'positions','horizontalalignment','center','fontsize',Fs-1,'color',sb{1},'verticalalignment','top')
    
    set(gca,'YAxisLocation','right')
    
    text(2,1.15,'\color[rgb]{0.7686 0.3059 0.3216}stimulus decoding','horizontalalignment','center','fontsize',9)
    text(2,0,sprintf('30 subsamples of\n%i%s stimulus neurons\\color[rgb]{0 0 0}',Nselunits,ss{4}),'horizontalalignment','center','fontsize',8)
    
    ps = arrayfun(@(i) signrank(decmat2(:,i),.25),1:3);
    label_ps(ps,1.03)
    
    arrayfun(@(x) set(x,'Color',sb{4}),B(:,1))
    arrayfun(@(x) set(x,'Color',sb{3}),B(:,2))
    arrayfun(@(x) set(x,'Color',sb{1}),B(:,3))
    
    text(3.45,.21,'chance','fontsize',8,'horizontalalignment','right')
    xlim([.5 3.5])
    set(gca,'Position',AX2_POS)
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 0.95,'pooled decoding during picture presentations','horizontalalignment','center','fontsize',9)
end