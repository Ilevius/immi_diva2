close all;
clc;
fieldName = fileread('C:\Users\tiama\OneDrive\??????? ????\IMMI\DIVA2\data\flat\1 upp\fieldname.txt');
waveLength = load('C:\Users\tiama\OneDrive\??????? ????\IMMI\DIVA2\data\flat\1 upp\waveLength.txt');
waveLength = waveLength(1,1);

fieldName = string(fieldName);
fieldName = strip(fieldName);
f1 = figure;
% f2 = figure;
% f3 = figure;
f1.WindowState = 'maximized';
% f2.WindowState = 'maximized';
% f3.WindowState = 'maximized';
%                               The    A B S O L U T E   value
figure(f1);
v = load('C:\Users\tiama\OneDrive\??????? ????\IMMI\DIVA2\data\flat\1 upp\integral_abs.txt');
dots = v(:, 1)/waveLength;
psi = v(1, 3);
psi = 'R/\lambda, \psi = '+string(psi)+'\circ';
graph1u = v(:, 2);

v = load('C:\Users\tiama\OneDrive\??????? ????\IMMI\DIVA2\data\flat\2 wpp\integral_abs.txt');
graph1w = v(:, 2);

graph1 = sqrt(graph1u.^2 + graph1w.^2);


v = load('C:\Users\tiama\OneDrive\??????? ????\IMMI\DIVA2\data\flat\1 upp\asymptotics_abs.txt');
graph2u = v(:, 2);
v = load('C:\Users\tiama\OneDrive\??????? ????\IMMI\DIVA2\data\flat\2 wpp\asymptotics_abs.txt');
graph2w = v(:, 2);

graph2 = sqrt(graph2u.^2 + graph2w.^2);

plot(dots, graph1, dots, graph2, 'x', 'LineWidth', 3, 'MarkerSize', 10);
ylim([0 inf]);
%xlim([10 inf]);
xlabel(psi, 'FontSize',36);
%ylabel('$|u|$, m', 'FontSize',36, 'Interpreter', 'Latex');
%legend(' |u_{-}^{P}| integral', '|A(\psi,R)| asymptotics');

legItem1 = '$|' + fieldName + '|$' + ' integral';
legItem2 = '$|' + fieldName + '|$' + ' asymptotics';
leg1 = legend(legItem1, legItem2);
set(leg1, 'Interpreter', 'Latex');


set(gca, 'FontSize',36);
grid on;
set(gcf,'color','w');
% 
% %                               The   R E A L   part
% figure(f2);
% v = load('C:\Users\tiama\OneDrive\??????? ????\IMMI\DIVA2\data\integral_real.txt');
% dots = v(:, 1);
% graph1 = v(:, 2);
% 
% v = load('C:\Users\tiama\OneDrive\??????? ????\IMMI\DIVA2\data\asymptotics_real.txt');
% graph2 = v(:, 2);
% 
% plot(dots, graph1, dots, graph2, 'x', 'LineWidth', 3, 'MarkerSize', 10);
% xlabel(psi, 'FontSize',36);
% %ylabel('$|u|$, m', 'FontSize',36, 'Interpreter', 'Latex');
% legend(' Re(u_{-}^{P}) integral', 'Re(A(\psi,R)) asymptotics');
% 
% set(gca, 'FontSize',36);
% grid on;
% set(gcf,'color','w');
% 
% %                             The   I M A G I N A R Y   part
% figure(f3);
% v = load('C:\Users\tiama\OneDrive\??????? ????\IMMI\DIVA2\data\integral_imag.txt');
% dots = v(:, 1);
% graph1 = v(:, 2);
% 
% v = load('C:\Users\tiama\OneDrive\??????? ????\IMMI\DIVA2\data\asymptotics_imag.txt');
% graph2 = v(:, 2);
% 
% plot(dots, graph1, dots, graph2, 'x', 'LineWidth', 3, 'MarkerSize', 10);
% xlabel(psi, 'FontSize',36);
% %ylabel('$|u|$, m', 'FontSize',36, 'Interpreter', 'Latex');
% legend(' Im(u_{-}^{P}) integral', 'Im(A(\psi,R)) asymptotics');
% 
% set(gca, 'FontSize',36);
% grid on;
% set(gcf,'color','w');
% 
% 
% 
% 
