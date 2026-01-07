%%
clear all
load fig3Fdata
cccc = get_seaborn;
f = figure('Position',[680   558   540   420]);

H1 = cplot(selected_resp_times,patwise_context_decoding,{'color',cccc{3}}); 
hold on,plot([-2000 1000],[.2 .2],'k--')
H2 = cplot(selected_resp_times,patwise_im1_decoding,{'color',cccc{8}});
H3 = cplot(selected_resp_times,patwise_im2_decoding,{'color',cccc{5}});
hold on,plot([-2000 1000],[.25 .25],'k--')


ylabel('decoding accuracy','FontSize',12)
xlabel('time until decision in ms')

AX=gca;

hold on 
% contextdec_stats: permttest stats of context decoding vs chance decoding (0.2)
for u = 1:numel(contextdec_stats.froms)
    H_significant = plot(selected_resp_times([(contextdec_stats.froms(u)) (contextdec_stats.tos(u))]),[1 1]*AX.YLim(1)*1+.025,'LineWidth',2,'Color',cccc{3});
end

% im1dec_stats: permttest stats of im1 decoding vs chance decoding (0.25)
for u = 1:numel(im1dec_stats.froms)
    H_significant = plot(selected_resp_times([(im1dec_stats.froms(u)) (im1dec_stats.tos(u))]),[1 1]*AX.YLim(1)*1.03+.025,'LineWidth',2,'Color',cccc{8});
end

% im2dec_stats: permttest stats of im2 decoding vs chance decoding (0.25)
for u = 1:numel(im2dec_stats.froms)
    H_significant = plot(selected_resp_times([(im2dec_stats.froms(u)) (im2dec_stats.tos(u))]),[1 1]*AX.YLim(1)*1.06+.025,'LineWidth',2,'Color',cccc{5});
end

plot([0 0],AX.YLim,'-','Linewidth',2,'color',cccc{4})
text(-470,AX.YLim(2)*.96,'decision','fontsize',12)

[hleg1, hobj1] = legend([H1.patch,H2.patch,H3.patch],{sprintf('context (CN, N = %i)',size(patwise_context_decoding,1)),sprintf('pic1 (SN, N = %i)',size(patwise_im1_decoding,1)),sprintf('pic2 (SN, N = %i)',size(patwise_im2_decoding,1))},'location','Northwest','fontsize',11);

XL = get(gca,'XLim');
YL = get(gca,'YLim');


plot([0 0],AX.YLim,'-','Linewidth',2,'color',cccc{4})
set(gca,'YAxisLocation','right')

%
text(250,.193,'chance (context)')
text(250,.193+.05,'chance (pictures)')

%% Code linear SVM context decoding ( LIBSVM library (v3.24) )
function [correct_preds] = context_dec_s_norm(x, y)

    if ~exist('libsetting', 'var') || isempty(libsetting)
        % Set libsetting for a linear decoder
        libsetting = '-s 0 -t 0 -c 1';
    end

    nRep = 5;
    nKs = 5;

    Nquestion = numel(unique(y));

    correct_preds = zeros(nRep, Nquestion);

    for rep = 1:nRep
        indices = crossvalind('Kfold', y, nKs);
        all_pred = nan(1, numel(y));

        for k = 1:nKs
            test = (indices == k);
            train = ~test;

            % Calculate mean and std from training data
            M = mean(x(train, :), 1);
            STD = std(x(train, :), [], 1);

            % Normalize training and test data
            x_train_norm = (x(train, :) - M) ./ STD;
            x_test_norm = (x(test, :) - M) ./ STD;

            % Train the SVM model on normalized data
            model = svmtrain(y(train)', x_train_norm, libsetting);

            % Predict using the normalized test data
            [predict_label, accuracy, dec_values] = svmpredict(y(test)', x_test_norm, model);

            all_pred(test) = predict_label';
        end

        % Calculate correct predictions for each group in y
        correct_preds(rep, :) = grpstats(all_pred == y, y)';
    end
end