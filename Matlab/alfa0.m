close all;
clc;

v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\alfa0.txt');
Rs = v(:, 1);
alfas = abs(v(:, 2));
kappa1 = v(:, 3);



v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\kappa2.txt');
kappa2 = v(:, 2);
psi = v(1, 1);
psi = 'R, \psi = '+string(psi)+'\circ';


plot(Rs, kappa1, Rs, kappa2, Rs, alfas, 'lineWidth', 3);
xlabel(psi, 'FontSize',36);
%ylabel('$|u|$, m', 'FontSize',36, 'Interpreter', 'Latex',);
legend(' \kappa_1', '\kappa_2', '\alpha_0');

set(gca, 'FontSize',36);
grid on;
set(gcf,'color','w');
