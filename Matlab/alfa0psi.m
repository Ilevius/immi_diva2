close all;
clc;

v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\alfa0psi.txt');
R = v(1, 1);
R = '\psi \circ, R = '+string(R)+'mm';
alfas = abs(v(:, 2));
kappa1 = v(:, 3);



v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\kappa2psi.txt');
kappa2 = v(:, 2);
psi = v(:, 1);


plot(psi, kappa1, psi, kappa2, psi, alfas, 'lineWidth', 3);
xlabel(R, 'FontSize',36);
%ylabel('$|u|$, m', 'FontSize',36, 'Interpreter', 'Latex',);
legend(' \kappa_1', '\kappa_2', '|\alpha_0|');

set(gca, 'FontSize',36);
grid on;
set(gcf,'color','w');
