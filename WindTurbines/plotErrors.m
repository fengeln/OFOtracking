clc
close all
clear all
colorvec = [
    254 97 0;
    120  94 240;
    220  38 127
] / 255;

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

eta = 0.1;
filename = sprintf("data/simulation_OFO_eta_%.2f.csv", eta);
OFOdata = readtable(filename);
u1      = OFOdata.u1;
u2      = OFOdata.u2;
uopt1   = OFOdata.uopt1;
uopt2   = OFOdata.uopt2;
windDir = OFOdata.windDir;
tsec    = OFOdata.tsec;
stepsize   = OFOdata.eta(1);
error1 = zeros(length(u1));
for tid = 1:length(u1)
    e = [u1(tid);u2(tid)] - [uopt1(tid);uopt2(tid)];
    error1(tid) =  sqrt(e' * D * e);
end
eta = 0.05;
filename = sprintf("data/simulation_OFO_eta_%.2f.csv", eta);
OFOdata = readtable(filename);
u1      = OFOdata.u1;
u2      = OFOdata.u2;
uopt1   = OFOdata.uopt1;
uopt2   = OFOdata.uopt2;
windDir = OFOdata.windDir;
tsec    = OFOdata.tsec;
stepsize   = OFOdata.eta(1);
error1 = zeros(length(u1),1);
for tid = 1:length(u1)
    e = [u1(tid);u2(tid)] - [uopt1(tid);uopt2(tid)];
    error1(tid) =  sqrt(e' * D * e);
end
eta = 0.155;
filename = sprintf("data/simulation_OFO_eta_%.2f.csv", eta);
OFOdata = readtable(filename);
u1      = OFOdata.u1;
u2      = OFOdata.u2;
uopt1   = OFOdata.uopt1;
uopt2   = OFOdata.uopt2;
windDir = OFOdata.windDir;
tsec    = OFOdata.tsec;
stepsize   = OFOdata.eta(1);
error2 = zeros(length(u1),1);
for tid = 1:length(u1)
    e = [u1(tid);u2(tid)] - [uopt1(tid);uopt2(tid)];
    error2(tid) =  sqrt(e' * D * e);
end
eta = 0.231;
filename = sprintf("data/simulation_OFO_eta_%.2f.csv", eta);
OFOdata = readtable(filename);
u1      = OFOdata.u1;
u2      = OFOdata.u2;
uopt1   = OFOdata.uopt1;
uopt2   = OFOdata.uopt2;
windDir = OFOdata.windDir;
tsec    = OFOdata.tsec;
stepsize   = OFOdata.eta(1);
error3 = zeros(length(u1),1);
for tid = 1:length(u1)
    e = [u1(tid);u2(tid)] - [uopt1(tid);uopt2(tid)];
    error3(tid) =  sqrt(e' * D * e);
end

figure(1)
plot(etavec,totalerror)
grid on

figure(2)
plot(tsec/60,error1,'Color',colorvec(3,:), 'LineWidth',2)
hold on
grid on
plot(tsec/60,error2,'Color',colorvec(1,:), 'LineWidth',2)
plot(tsec/60,error3,'Color',colorvec(2,:), 'LineWidth',2)
ylabel('$\|\epsilon_k\|_D$','Interpreter','latex','FontSize',12)
xlabel('time in minutes','Interpreter','latex','FontSize',12)
legend({'$\eta = 0.05$','$\eta_{\mathrm{opt}} = 0.155$','$\eta^\ast = 0.231$'},'Interpreter', 'latex','FontSize', 12,'Location', 'northwest');
fig = gcf;
fig.Position = [100 100 400 200];

% -------------------------------------------------------------------------
%% Export plots
% -------------------------------------------------------------------------
figure(2)
exportgraphics(gcf,'figures/WindOFOErrors.pdf', ContentType='vector')