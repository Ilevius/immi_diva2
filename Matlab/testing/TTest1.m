close all;
clc;
v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\TTestGr1.txt');
alfas = v(:, 1);
left_part = v(:, 2);
right_part = v(:,3);

plot(alfas, left_part, alfas, right_part, 'x', 'LineWidth', 3, 'MarkerSize', 10);
legend('Re(\mu_1(U_0^{\prime}(\alpha, -h)+U_-^{\prime}(\alpha, -h)-i\alpha(W_0(\alpha, -h)+W_-(\alpha, -h))))', 'Re(\mu_2(U_+^{\prime}(\alpha, -h)-i\alpha W_+(\alpha, -h)))');

set(gca, 'FontSize',36);
xlabel('\alpha', 'FontSize',36);
%ylabel('$|u|$, m', 'FontSize',36, 'Interpreter', 'Latex');

grid on;
set(gcf,'color','w');