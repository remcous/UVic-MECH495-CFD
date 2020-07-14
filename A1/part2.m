% Initialize constants from problem definition
W = 0.5;
D = 1;
Bi = 1;
alpha = [0.1 0.2 0.3];
epsilon = 1e-5;

% Set function handles
t=@(xi) W/2*(1-xi)^2;
Ac = @(xi_j) D*W*(1-xi_j)^2;
d_Ac =@(xi_j) -2*D*W*(1-xi_j);
d_As1 =@(xi_j, d_xi) 2*D/d_xi * sqrt((t(xi_j)-t(xi_j+d_xi))^2+d_xi^2);
d_As2 =@(xi_j, d_xi) W*d_xi*(-4+4*xi_j+2*d_xi);
d_As =@(xi_j, d_xi) d_As1(xi_j, d_xi) + d_As2(xi_j, d_xi);
p =@(xi_j) d_Ac(xi_j)/Ac(xi_j);
r =@(xi_j, d_xi, alpha, y_old) -1/((1-alpha * y_old)*Ac(xi_j))*d_As(xi_j, d_xi);
a =@(xi_j, d_xi, alpha, y_old) -2+d_xi^2*r(xi_j,d_xi, alpha, y_old);
r_linear =@(xi_j, d_xi) -Bi/Ac(xi_j)*d_As(xi_j, d_xi);
a_linear =@(xi_j, d_xi) -2+d_xi^2*r_linear(xi_j,d_xi);
b =@(xi_j, d_xi) 1+d_xi/2*p(xi_j);
c =@(xi_j, d_xi) 1-d_xi/2*p(xi_j);

% run linear solver for 50 and 400 mesh intervals to get initial y value
[xi50, y50_initial] = solver1D(a_linear,b,c,Bi,50);
[y_50_1, e50_1] = solver1Dnonlinear(a, b, c, alpha(1), xi50, y50_initial, 50, epsilon);
[y_50_2, e50_2] = solver1Dnonlinear(a, b, c, alpha(2), xi50, y50_initial, 50, epsilon);
[y_50_3, e50_3] = solver1Dnonlinear(a, b, c, alpha(3), xi50, y50_initial, 50, epsilon);
[xi400, y400_initial] = solver1D(a_linear,b,c,Bi,400);
[y_400_1, e400_1] = solver1Dnonlinear(a, b, c, alpha(1), xi400, y400_initial, 400, epsilon);
[y_400_2, e400_2] = solver1Dnonlinear(a, b, c, alpha(2), xi400, y400_initial, 400, epsilon);
[y_400_3, e400_3] = solver1Dnonlinear(a, b, c, alpha(3), xi400, y400_initial, 400, epsilon);

% Top plot - full view
subplot(3,1,1)
hold on;    %only use one plot window
p1=plot(xi50,y_50_1,'--r','LineWidth',3);
p2=plot(xi400,y_400_1,'-b','LineWidth',2);
legend('50 Intervals', '400 Intervals');
xlabel('\xi'), ylabel('y');
title("Non-linear Thermal Conductivitiy with \alpha=0.1");
hold off;

subplot(3,1,2)
hold on;    %only use one plot window
plot(xi50,y_50_2,'--r','LineWidth',3);
plot(xi400,y_400_2,'-b','LineWidth',2);
legend('50 Intervals', '400 Intervals');
xlabel('\xi'), ylabel('y');
title("Non-linear Thermal Conductivitiy with \alpha=0.2");
hold off;

subplot(3,1,3)
hold on;    %only use one plot window
plot(xi50,y_50_3,'--r','LineWidth',3);
plot(xi400,y_400_3,'-b','LineWidth',2);
legend('50 Intervals', '400 Intervals');
xlabel('\xi'), ylabel('y');
title("Non-linear Thermal Conductivitiy with \alpha=0.3");
hold off;

% saves figure
saveas(gcf, 'question2_iteration.png');

figure;     % creates new figure
plot(xi400,y_400_1,'-r','LineWidth',2);
hold on;    %only use one plot window
plot(xi400,y_400_2,'-g','LineWidth',2);
plot(xi400,y_400_3,'-b','LineWidth',2);
legend('\alpha=0.1', '\alpha=0.2', '\alpha=0.3');
xlabel('\xi'), ylabel('y');
title("Non-linear Thermal Conductivitiy with varying \alpha");
hold off;

% saves figure
saveas(gcf, 'question2_alpha.png');