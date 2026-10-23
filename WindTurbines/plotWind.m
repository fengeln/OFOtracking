clc
close all 
clear all
addpath("..")
x = readmatrix('data/wind_1s.csv');
blueRGB = [100 143 255]/255;

t=1:1:length(x);
figure();
plot(t./60,x,'-', 'LineWidth', 2, 'Color', blueRGB)
grid on
xlabel('time in minutes','Interpreter', 'latex','FontSize', 12)
ylabel('$\theta$ in degrees', 'Interpreter', 'latex','FontSize', 12)
xlim([0 100])
ylim([254 286])
fig = gcf;
fig.Position = [100 100 400 120];
exportgraphics(gcf,'figures/windProfile.pdf', ContentType='vector')