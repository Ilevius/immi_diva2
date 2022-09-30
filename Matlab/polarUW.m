close all;
clc;

fieldName = fileread('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\5 usp\fieldname.txt');
fieldName = string(fieldName);
fieldName = strip(fieldName);

v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\5 usp\integral_abs.txt');
psis = v(:, 3);
v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\5 usp\integral_abs.txt');
R = v(1, 1);
R = string(R);



v1 = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\5 usp\integral_abs.txt');
v2 = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\6 wsp\integral_abs.txt');
graph1 = sqrt(v1(:, 2).^2 + v2(:, 2).^2);

v1 = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\5 usp\asymptotics_abs.txt');
v2 = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\6 wsp\asymptotics_abs.txt');
graph2 = sqrt(v1(:, 2).^2 + v2(:, 2).^2);

psis = psis/180*pi;

polarplot(psis, graph1, psis, graph2, 'x', 'LineWidth', 3, 'MarkerSize', 10);
thetalim([0 180]);


% xlabel('\psi, R  =  '+R+' mm', 'FontSize',36);
% %ylabel('$|u|$, m', 'FontSize',36, 'Interpreter', 'Latex');


legItem1 = '$|' + fieldName + '|$' + ' integral';
legItem2 = '$|' + fieldName + '|$' + ' asymptotics';
leg1 = legend(legItem1, legItem2);
set(leg1, 'Interpreter', 'Latex');

set(gca, 'FontSize',32);
grid on;
set(gcf,'color','w');
