close all;
clc;
v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\TTestGr2.txt');
alfas = v(:, 1);
left_part = v(:, 2);
right_part = v(:,3);

plot(alfas, left_part, alfas, right_part, 'x', 'LineWidth', 3, 'MarkerSize', 10);
%legend('Re(-i\alpha\lambda_1 (U_0(\alpha, -h)+U_-(\alpha, -h)) + (\lambda_1+2\mu_1)(W_0^{\prime}(\alpha, -h)+W_-^{\prime}(\alpha, -h)))', 'Re(-i\alpha\lambda_2 U_+(\alpha, -h) + (\lambda_2+2\mu_2)W_+^{\prime}(\alpha, -h))');
leg1 = legend('$|T_z(V_0(\alpha, -h) + V_-^1(\alpha, -h))|$', '$|T_z(V_+^1(\alpha, -h))|$');
set(leg1, 'Interpreter', 'Latex');


set(gca, 'FontSize',36);
xlabel('\alpha', 'FontSize',36);
%ylabel('$|u|$, m', 'FontSize',36, 'Interpreter', 'Latex');

grid on;
set(gcf,'color','w');