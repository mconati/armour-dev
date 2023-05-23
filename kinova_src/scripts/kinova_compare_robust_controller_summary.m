close all; clear; clc;

model_uncertainty = [0, 3, 5, 10, 25, 50];

x = 1:length(model_uncertainty);
middle = zeros(7, 2, length(model_uncertainty));
errhigh = zeros(7, 2, length(model_uncertainty));
errlow = zeros(7, 2, length(model_uncertainty));  

for i = 1:length(model_uncertainty)
    data = load(['robust_input_', num2str(model_uncertainty(i)), '.mat']);
    middle(:,1,i) = median(abs(data.robust_out_armour_final), 2);
    errhigh(:,1,i) = max(abs(data.robust_out_armour_final), [], 2) - middle(:,1,i);
    errlow(:,1,i) = middle(:,1,i) - min(abs(data.robust_out_armour_final), [], 2);
    middle(:,2,i) = median(abs(data.robust_out_althoff_final), 2);
    errhigh(:,2,i) = max(abs(data.robust_out_althoff_final), [], 2) - middle(:,2,i);
    errlow(:,2,i) = middle(:,2,i) - min(abs(data.robust_out_althoff_final), [], 2);
end

figure;
for i = 1:7
    subplot(3,3,i);
    hBar = bar(x,squeeze(middle(i,:,:)));   

    hold on
%     er = errorbar(x,squeeze(middle(i,:,:)),squeeze(errlow(i,:,:)),squeeze(errhigh(i,:,:)));    
% %     er.Color = [0 0 0];                            
%     er(1).LineStyle = 'none';
%     er(2).LineStyle = 'none';
    for j = 1:size(squeeze(middle(i,:,:)),1)
        ctr(j,:) = bsxfun(@plus, hBar(j).XData, hBar(j).XOffset');   % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
        ydt(j,:) = hBar(j).YData;                                     % Individual Bar Heights
    end
    hold on
    er = errorbar(ctr(1,:), ydt(1,:), squeeze(errlow(i,1,:)), squeeze(errhigh(i,1,:)));
    er.LineStyle = 'none';
    er.Color = [0, 0.4470, 0.7410];
    er = errorbar(ctr(2,:), ydt(2,:), squeeze(errlow(i,2,:)), squeeze(errhigh(i,2,:)));
    er.LineStyle = 'none';
    er.Color = [0.8500, 0.3250, 0.0980];

    ylabel('robust input (N*m)');
    xlabel('model uncertainty');
    set(gca,'XTickLabel',{'0%','3%','5%','10%','25%','50%'});

    if i == 1
        legend('armour', 'althoff');
    end
end

sgtitle('robust input comparison (median + min/max error bar)');