close all;
clc;
v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\dots.txt');
x = v(:, 1);
z = v(:, 2);

v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\ReDots.txt');
Rex = v(:, 1);
Rez = v(:, 2);

v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\ReDots1.txt');
Re1x = v(:, 1);
Re1z = v(:, 2);


plot(x, z, Rex, Rez, 'x', Re1x, Re1z, 'o', 'LineWidth', 3, 'MarkerSize', 10);
legend(' old dots', 'new dots from 2h', 'dots from h');

set(gca, 'FontSize',36);
xlabel('x', 'FontSize',36);
ylabel('$z$, mm', 'FontSize',36, 'Interpreter', 'Latex');

grid on;
set(gcf,'color','w');