x = [26*51 51*101 201*101];
jacobi = [684 1720 4285];
gauss = [558 1418 3577];
sor = [424 1096 2770];
adi = [243 685 1722];

hold on;
plot(x, jacobi, x, gauss, x, sor, x, adi);
legend('Jacobi', 'Gauss-Seidel', 'SOR', 'ADI', 'Location', 'northwest');
xlabel('Number of Mesh Nodes'), ylabel('Number of Iterations');
set(gca, 'XScale', 'log');
title("Number of Iterations to Converge by Mesh Size");
grid on;
saveas(gcf,'convergencePlot.png'); close(gcf);