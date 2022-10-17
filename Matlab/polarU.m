close all;
clc;

fieldName = fileread('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\fieldname.txt');
fieldName = string(fieldName);
fieldName = strip(fieldName);

v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\integral_abs.txt');
psis = v(:, 2);
v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\integral_abs.txt');
R = v(1, 1);
R = string(R);



v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\integral_abs.txt');
graph1 = v(:, 3);

v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\asymptotics_abs.txt');
graph2 = v(:, 3);

%psis = psis/180*pi;

polarplot(psis, graph1, psis, graph2, '--', 'LineWidth', 3, 'MarkerSize', 10);
%polarplot(psis, graph1, 'LineWidth', 3, 'MarkerSize', 10);
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
