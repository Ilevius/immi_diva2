close all;
clc;

fieldName = fileread('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\polar\3 ups\fieldname.txt');
fieldName = string(fieldName);
fieldName = strip(fieldName);

v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\polar\3 ups\integral_abs.txt');
psis = v(:, 3);
v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\polar\3 ups\integral_abs.txt');
R = v(1, 1);
R = string(R);



v1 = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\polar\3 ups\integral_abs.txt');
v2 = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\polar\4 wps\integral_abs.txt');
graph1 = sqrt(v1(:, 2).^2 + v2(:, 2).^2);

v1 = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\polar\3 ups\asymptotics_abs.txt');
v2 = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\polar\4 wps\asymptotics_abs.txt');
graph2 = sqrt(v1(:, 2).^2 + v2(:, 2).^2);

psis = psis/180*pi;

m1 = max(graph1);
m2 = max(graph2);
maxi = max(m1, m1)

graph1 = graph1/maxi;
graph2 = graph2/maxi;

polarplot(psis, graph1, psis, graph2, '--', 'LineWidth', 1, 'MarkerSize', 28);
thetalim([0 180]);


% xlabel('\psi, R  =  '+R+' mm', 'FontSize',36);
% %ylabel('$|u|$, m', 'FontSize',36, 'Interpreter', 'Latex');


% legItem1 = '$|' + fieldName + '|$' + ' integral';
% legItem2 = '$|' + fieldName + '|$' + ' asymptotics';
% leg1 = legend(legItem1, legItem2);
% set(leg1, 'Interpreter', 'Latex');
grid on;
set(gca, 'FontSize',32);
set(gca, 'GridAlpha', 1);

set(gcf,'color','w');
