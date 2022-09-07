close all;
clc;

v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\integral_abs.txt');
psis = v(:, 3);
v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\integral_abs.txt');
R = v(1, 1);
R = string(R);


f1 = figure;
f2 = figure;
f3 = figure;
f1.WindowState = 'maximized';
f2.WindowState = 'maximized';
f3.WindowState = 'maximized';
%                               The    A B S O L U T E   value
figure(f1);

v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\integral_abs.txt');
graph1 = v(:, 2);

v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\asymptotics_abs.txt');
graph2 = v(:, 2);

plot(psis, graph1, psis, graph2, 'x', 'LineWidth', 3, 'MarkerSize', 10);
xlabel('\psi, R  =  '+R+' mm', 'FontSize',36);
%ylabel('$|u|$, m', 'FontSize',36, 'Interpreter', 'Latex');
legend(' |u_{-}^{P}| integral', '|A(\psi,R)| asymptotics');

set(gca, 'FontSize',36);
grid on;
set(gcf,'color','w');

%                               The   R E A L   part
figure(f2);

v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\integral_real.txt');
graph1 = v(:, 2);

v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\asymptotics_real.txt');
graph2 = v(:, 2);

plot(psis, graph1, psis, graph2, 'x', 'LineWidth', 3, 'MarkerSize', 10);
xlabel('\psi, R  =  '+R+' mm', 'FontSize',36); 
%ylabel('$|u|$, m', 'FontSize',36, 'Interpreter', 'Latex');
legend(' Re(u_{-}^{P}) integral', 'Re(A(\psi,R)) asymptotics');

set(gca, 'FontSize',36);
grid on;
set(gcf,'color','w');

%                             The   I M A G I N A R Y   part
figure(f3);

v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\integral_imag.txt');
graph1 = v(:, 2);

v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\asymptotics_imag.txt');
graph2 = v(:, 2);

plot(psis, graph1, psis, graph2, 'x', 'LineWidth', 3, 'MarkerSize', 10);
xlabel('\psi, R  =  '+R+' mm', 'FontSize',36);
%ylabel('$|u|$, m', 'FontSize',36, 'Interpreter', 'Latex');
legend(' Im(u_{-}^{P}) integral', 'Im(A(\psi,R)) asymptotics');

set(gca, 'FontSize',36);
grid on;
set(gcf,'color','w');




