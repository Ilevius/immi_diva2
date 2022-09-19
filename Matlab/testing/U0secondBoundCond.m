close all;
clc;
v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\U0testGr2.txt');
alfas = v(:, 1);
left_part = v(:, 2);
right_part = v(:,3);
title('Plots must coinside')
plot(alfas, left_part, alfas, right_part, 'x', 'LineWidth', 3, 'MarkerSize', 10);
%legend('Re(\lambda_1(-i\alpha U_0)+ (\lambda_1+2\mu_1) W_0^{\prime})', 'Q_2 (= 1)');

leg1 = legend('$|F_x[\sigma_z]|$', '$|Q_2| = 1$');
set(leg1, 'Interpreter', 'Latex');

set(gca, 'FontSize',36);
xlabel('\alpha, z=0', 'FontSize',36);
%ylabel('$|u|$, m', 'FontSize',36, 'Interpreter', 'Latex');

grid on;
set(gcf,'color','w');