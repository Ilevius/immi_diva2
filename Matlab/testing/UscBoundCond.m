close all;
clc;
v = load('C:\Users\tiama\OneDrive\??????? ????\IMMI\DIVA2\data\UscTest.txt');
alfas = v(:, 1);
left_part = v(:, 2);
right_part = v(:,3);

plot(alfas, left_part, alfas, right_part, 'x', 'LineWidth', 3, 'MarkerSize', 10);

leg1 = legend('$|U_0(\alpha, -h) + U_-^1(\alpha, -h)|$', '$|U_+^1(\alpha, -h)|$');
set(leg1, 'Interpreter', 'Latex');

set(gca, 'FontSize',36);
xlabel('\alpha, z=-h', 'FontSize',36);
%ylabel('$|u|$, m', 'FontSize',36, 'Interpreter', 'Latex');

grid on;
set(gcf,'color','w');