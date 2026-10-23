clear; close all; clc;
colorvec = [254 97 0;120  94 240;220  38 127]/ 255;
% Eigenvalue bounds
lambda_min = 1;
lambda_max = 2;
% time-varying parameter
N = 1000;
t = 0:1:N;
theta1 = (0.001*sin(0.99*pi*t) + 0.5*cos(0.01*pi*t))*1.28; 
theta2 = 0.01*sin(0.9*pi*t) + 0.1*cos(0.01*pi*t);
[J1,deltamax,etavec,etaopt1, Jopt1] = calcJ(N,theta1,lambda_min,lambda_max);
[J2,deltamax2,etavec2,etaopt2, Jopt2] = calcJ(N,theta2,lambda_min,lambda_max);

%% Tracking error bound r(eta)
%-------------------------------------------------------------
r = zeros(size(etavec));
for i = 1:length(etavec)
    rho = max(abs(1 - etavec(i)*lambda_max),abs(1 - etavec(i)*lambda_min));
    r(i) = deltamax/(1 - rho);
end
eta_opt = 2/(lambda_min + lambda_max);
rho_star = max(abs(1 - eta_opt*lambda_max),abs(1 - eta_opt*lambda_min));
r_star = deltamax/(1 - rho_star);
%-------------------------------------------------------------
%% Plot
load('data/error1.mat')
load('data/error2.mat')
etavecshort = 0.01:0.01:1.0;
figure(1);
hold on
box on
grid on
% Main curve
plot(etavec, r,'LineWidth', 2, 'Color', colorvec(1,:));
plot(etavec, J1./N, 'LineWidth', 2, 'Color', colorvec(2,:));
plot(etavec, J2./N, 'LineWidth', 2, 'Color', colorvec(3,:));

plot(etavecshort,totalerror2,'--','LineWidth', 2, 'Color', colorvec(2,:))
plot(etavecshort,totalerror1,'--','LineWidth', 2, 'Color', colorvec(3,:))
% Vertical dashed line
%xline(eta_opt, '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
%xline(etaopt1, '--', 'LineWidth', 2,'Color', [0.5 0.5 0.5]);
%xline(etaopt2, '--', 'LineWidth', 2,'Color', [0.5 0.5 0.5]);

% Marker at optimum
plot(eta_opt, r_star, 'o','MarkerSize',8,'MarkerFaceColor', colorvec(1,:),'MarkerEdgeColor', 'k');
plot(etaopt1, Jopt1/N, 'o','MarkerSize',8,'MarkerFaceColor',colorvec(2,:),'MarkerEdgeColor', 'k');
plot(etaopt2,Jopt2/N,'o','MarkerSize',8,'MarkerFaceColor',colorvec(3,:),'MarkerEdgeColor', 'k');
% Labels
xlabel('stepsize $\eta$','Interpreter','latex','FontSize', 15);

ylabel('(-)','Interpreter','latex','FontSize',15);

% Legend
legend({'$r(\eta)$','$\hat{J}(\eta)$','$\tilde{J}(\eta)$','$\sum\|\hat{\epsilon}_k\|_D$','$\sum\|\tilde{\epsilon}_k\|_D$'},'Interpreter', 'latex','FontSize',13,'Location','northwest');

% Axes styling
set(gca,'FontSize', 13,'LineWidth', 1.2);
xlim([min(etavec), max(etavec)])
ylim([0, 15])
xlim([0.0, 0.99])
fig = gcf;
fig.Position = [100 100 600 250];
%exportgraphics(gcf,'figures/ExampleEtaOpt.pdf', ContentType='vector')

function [J, deltamax, etavec,etaopt,Jopt] = calcJ(N,theta,lambda_min,lambda_max)
    omegavec = linspace(-pi,pi,N);
    ustar = [theta;2*theta];
    delta = ustar(:,2:end)-ustar(:,1:end-1);
    deltamax = max(vecnorm(delta,2,1));
    Delta = fftshift(fft(delta,N,2),2);
    
    etavec = linspace(1e-4, 2/lambda_max - 1e-4, 1000);
    lambdavec = linspace(lambda_min, lambda_max, 1000);
    
    J = zeros(size(etavec));
    
    for etaid = 1:length(etavec)
        etaid
        eta = etavec(etaid);
        J_eta = 0;
        for omegaid = 1:length(omegavec)
            omega = omegavec(omegaid);
            denomvec = zeros(size(lambdavec));
            for lambdaid = 1:length(lambdavec)
                lambda = lambdavec(lambdaid);
                a = 1 - eta * lambda;
                denomvec(lambdaid) = abs(exp(1i*omega) - a)^2;
            end
            J_eta = J_eta + norm(Delta(:,omegaid))^2 / min(denomvec);
        end
        J(etaid) = J_eta;
    end
    
    % Find minimizing eta
    [Jopt, idx] = min(J);
    etaopt = etavec(idx);
end