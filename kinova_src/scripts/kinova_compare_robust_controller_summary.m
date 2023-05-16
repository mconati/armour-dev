close all; clear; clc;

model_uncertainty = [0, 2.5, 5, 10, 25, 50];

mean_robust_out_armour = zeros(7, length(model_uncertainty));
mean_robust_out_althoff = zeros(7, length(model_uncertainty));

for i = 1:length(model_uncertainty)
    data = load(['robust_input_', num2str(model_uncertainty(i)), '.mat']);
    mean_robust_out_armour(:,i) = data.mean_robust_out_armour_final;
    mean_robust_out_althoff(:,i) = data.mean_robust_out_althoff_final;
end

for i = 1:7
    subplot(3,3,i); hold on;
    plot(model_uncertainty, mean_robust_out_althoff(i,:), 'r*-');
    plot(model_uncertainty, mean_robust_out_armour(i,:), 'bo-');
    if i == 1
        legend('Althoff','Armour');
    end
    ylim([0, max([mean_robust_out_althoff(i,:), mean_robust_out_armour(i,:)])* 1.1]);
    xticks(model_uncertainty);
    xlabel('model uncertainty %');
    ylabel('mean robust input (N*m)');
    title(['joint ', num2str(i)]);
end
sgtitle('robust input comparison');