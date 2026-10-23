clear; close all; clc;
colorvec = [
    254 97 0;
    120  94 240;
    220  38 127
] / 255;
blueRGB = [100 143 255]/255;

% Eigenvalue bounds
lambda_min = 4.38;
lambda_max = 8.48;
%% Load wind data
wind_s = readmatrix('data/wind_1s.csv');
wind_s = wind_s(:);              % ensure column vector
wind_min = wind_s(1:60:60*100);
delta = zeros(length(wind_min)-1,1);
for i = 2:length(wind_min)
    delta(i-1) = wind_min(i) - wind_min(i-1);
end
deltamax = max(delta);
N = length(delta);

%% Frequency grid for DTFT
Ts = 1;   % sampling time
Nw = 2048;
w = linspace(-pi, pi, Nw);   % normalized rad/sample

%% Compute DTFT
Delta = zeros(1, Nw);
n = 0:N-1;

for k = 1:Nw
    Delta(k) = sum(delta.' .* exp(-1j * w(k) * n));
end

%% Convert to physical frequency (rad/sec)
Omega = w / Ts;

% eta range
eta_vals = linspace(1e-4, 2/lambda_max - 1e-4, 1000);
weights = abs(Delta);

% lambda grid for maximization
lambda_grid = linspace(lambda_min, lambda_max, 2000);
J = zeros(size(eta_vals));

for k = 1:length(eta_vals)
    k
    eta = eta_vals(k);
    J_eta = 0;
    for i = 1:length(Omega)
        w = Omega(i);
        denom_vals = zeros(size(lambda_grid));% compute max over lambda (via grid search)
        for j = 1:length(lambda_grid)
            lambda = lambda_grid(j);
            a = 1 - eta * lambda;
            denom_vals(j) = abs(exp(1i*w) - a)^2;
        end
        M_eta = 1 / min(denom_vals);
        J_eta = J_eta + weights(i)^2 * M_eta;
    end
    J(k) = J_eta;
end
% Find minimizing eta
[Jmin, idx] = min(J);
eta_star = eta_vals(idx);


% Tracking error bound r(eta)
r = zeros(size(eta_vals));
for i = 1:length(eta_vals)
    rho = max(abs(1 - eta_vals(i)*lambda_max), ...
              abs(1 - eta_vals(i)*lambda_min));
    r(i) = deltamax/(1 - rho);
end
eta_opt = 2 / (lambda_max + lambda_min);

% -------------------------------------------------------------------------
%% Calculate total error from OFO data
% -------------------------------------------------------------------------
etavec = 0.0:0.01:0.3;
D = [10 0;0 0.0008];
totalerror = zeros(length(etavec),1);

for etaid = 1:length(etavec)
    eta = etavec(etaid);
    filename = sprintf("data/simulation_OFO_eta_%.2f.csv", eta);
    OFOdata = readtable(filename);
    u1      = OFOdata.u1;
    u2      = OFOdata.u2;
    uopt1   = OFOdata.uopt1;
    uopt2   = OFOdata.uopt2;
    windDir = OFOdata.windDir;
    tsec    = OFOdata.tsec;
    stepsize   = OFOdata.eta(1);
    for tid = 1:length(u1)
        e = [u1(tid);u2(tid)] - [uopt1(tid);uopt2(tid)];
        totalerror(etaid) = totalerror(etaid) + sqrt(e' * D * e);
    end
end
% -------------------------------------------------------------------------
%% Plot results
% -------------------------------------------------------------------------
%% Plot spectrum
figure(1)
plot(Omega, abs(Delta), 'LineWidth', 1.5, 'Color', blueRGB);
grid on;
xlabel('frequency $\omega$', 'Interpreter','latex','FontSize', 12);
ylabel('spectrum', 'Interpreter','latex','FontSize', 12);
xlim([-pi pi])
fig = gcf;
fig.Position = [100 100 400 120];

%% Plot radius and J over stepsize
figure(2);
grid on;
hold on;
box on
plot(eta_vals, r,'LineWidth', 2, 'Color', colorvec(1,:));
plot(eta_vals, 0.00001*J, 'LineWidth', 2, 'Color', colorvec(2,:));
plot(etavec,totalerror, '--', 'LineWidth', 2, 'Color', colorvec(2,:))
% Vertical dashed line
%xline(eta_opt, '--', 'LineWidth', 2, 'Color', colorvec(1,:));
%xline(eta_star, '--', 'LineWidth', 2,'Color', colorvec(2,:));
% Marker at optimum
rho_star = max(abs(1 - eta_opt*lambda_max), ...
               abs(1 - eta_opt*lambda_min));
r_star = deltamax/(1 - rho_star);

plot(eta_opt, r_star, 'o', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', colorvec(1,:), ...
    'MarkerEdgeColor', 'k');
plot(eta_star, 0.00001*Jmin, 'o', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', colorvec(2,:), ...
    'MarkerEdgeColor', 'k');
% Legend
legend({'$r(\eta)$','$10^{-5}\cdot J(\eta)$','$\sum\|\epsilon_k\|_D$'},'Interpreter', 'latex','FontSize', 12,'Location', 'northeast');
ylim([0 30])
xlim([min(eta_vals), max(eta_vals)])
xlabel('stepsize $\eta$','Interpreter', 'latex','FontSize', 12);
ylabel('(-)','Interpreter', 'latex','FontSize', 12);
%set(gca, 'YScale', 'log')
fig = gcf;
fig.Position = [100 100 400 200];
% -------------------------------------------------------------------------
%% Export plots
% -------------------------------------------------------------------------
figure(1)
exportgraphics(gcf,'figures/WindSpectrum.pdf', ContentType='vector')

figure(2)
exportgraphics(gcf,'figures/WindStepsize.pdf', ContentType='vector')