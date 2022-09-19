close all;
clc;
v = load('C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\VscTest.txt');
alfas = v(:, 1);
left_part = v(:, 2);
right_part = v(:,3);

plot(alfas, left_part, alfas, right_part, 'x', 'LineWidth', 3, 'MarkerSize', 10);
%legend('Im(V_0(\alpha, -h) + V_-(\alpha, -h))', 'Im(V_+(\alpha, -h))');

leg1 = legend('$|V_0(\alpha, -h) + V_-^1(\alpha, -h)|$', '$|V_+^1(\alpha, -h)|$');
set(leg1, 'Interpreter', 'Latex');

set(gca, 'FontSize',36);
xlabel('\alpha, z=-h', 'FontSize',36);
%ylabel('$|u|$, m', 'FontSize',36, 'Interpreter', 'Latex');

grid on;
set(gcf,'color','w');