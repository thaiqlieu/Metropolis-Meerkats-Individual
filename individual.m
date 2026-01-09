close all;

N       = 100000;
burnin  = 10000;
sigmas   = [0.01,0.10,0.50];

acceptance_rates = zeros(length(sigmas),1);
samples_all      = cell(length(sigmas),1);


target = @(z) exp(-5*(z-0.5).^2) .* cos(3*pi*z) + 1;

for s = 1:length(sigmas)

    sigma = sigmas(s);
    accept = 0;
    x = zeros(N,1);
    x(1) = 0.5;
    
    for n = 2:N
        proposal = x(n-1) + sigma*randn;

        if proposal < 0 || proposal > 1
            x(n) = x(n-1);
            continue
        end
        alpha = min(1, target(proposal) / target(x(n-1)));

        if rand < alpha
            x(n) = proposal;
            accept = accept + 1;
        else
            x(n) = x(n-1);
        end
    end

    acceptance_rates(s) = accept / (N - 1);
    samples_all{s} = x(burnin+1:end);
    
    fprintf('Sigma = %.2f, Acceptance Rate = %.2f%%\n', sigma, 100*acceptance_rates(s));
end

% Figure 1: Acceptance Rates vs Proposal States
sigmas = sigmas(:);
T = table(sigmas, acceptance_rates,...
    'VariableNames', {'Sigma', 'AcceptanceRate_percent'});
figure;
uitable('Data', T{:,:},...
        'ColumnName', T.Properties.VariableNames,...
        'RowName', [],...
        'Units','Normalized',...
        'Position',[0 0 1 1]);
title('Acceptance Rates for Different Proposal Scales');

% Figure 2: Trace Plots
for i = 1:length(sigmas)
    figure;
    plot(1:5000, samples_all{i}(1:5000));
    xlabel('Iteration'); ylabel('X_n');
    title(['Trace Plot for \sigma = ', num2str(sigmas(i))]);
end

% Figure 3: Histograms vs Target Distribution
xx = linspace(0,1,400);
yy = target(xx);
yy = yy / trapz(xx, yy);

for i = 1:length(sigmas)
    figure;
    histogram(samples_all{i}, 40, 'Normalization','pdf');
    hold on;
    plot(xx, yy, 'k', 'LineWidth', 1.5);
    title(['\sigma = ', num2str(sigmas(i))]);
end
