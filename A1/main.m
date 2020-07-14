% Initialize constants from problem definition
W = 0.5;
D = 1;
Bi = 1;

% Set function handles
t=@(xi) W/2*(1-xi)^2;
Ac = @(xi_j) D*W*(1-xi_j)^2;
d_Ac =@(xi_j) -2*D*W*(1-xi_j);
d_As1 =@(xi_j, d_xi) 2*D/d_xi * sqrt((t(xi_j)-t(xi_j+d_xi))^2+d_xi^2);
d_As2 =@(xi_j, d_xi) W*d_xi*(-4+4*xi_j+2*d_xi);
d_As =@(xi_j, d_xi) d_As1(xi_j, d_xi) + d_As2(xi_j, d_xi);
p =@(xi_j) d_Ac(xi_j)/Ac(xi_j);
r =@(xi_j, d_xi) -Bi/Ac(xi_j)*d_As(xi_j, d_xi);
a =@(xi_j, d_xi) -2+d_xi^2*r(xi_j,d_xi);
b =@(xi_j, d_xi) 1+d_xi/2*p(xi_j);
c =@(xi_j, d_xi) 1-d_xi/2*p(xi_j);

% run solver for 50, 100, 200, and 400 mesh intervals
[xi50, y50] = solver1D(a,b,c,Bi,50);
[xi100, y100] = solver1D(a,b,c,Bi,100);
[xi200, y200] = solver1D(a,b,c,Bi,200);
[xi400, y400] = solver1D(a,b,c,Bi,400);

% Top plot - full view
subplot(2,1,1)
hold on;    %only use one plot window
plot(xi50,y50,'-r','LineWidth',5);
plot(xi100,y100,'--m','LineWidth',4);
plot(xi200,y200,'-.b','LineWidth',3);
plot(xi400,y400,':g','LineWidth',2);
legend('50 intervals', '100 intervals', '200 intervals', '400 intervals');
xlabel('\xi'), ylabel('y');
title("Non-dimensional temperature difference by mesh size");
hold off;

% Bottom plot - zoomed on convective tip
subplot(2,1,2);
hold on;
axis([0.95 1 0 0.01]);
plot(xi50,y50,'-r','LineWidth',5);
plot(xi100,y100,'--m','LineWidth',4);
plot(xi200,y200,'-.b','LineWidth',3);
plot(xi400,y400,':g','LineWidth',2);
legend('50 intervals', '100 intervals', '200 intervals', '400 intervals');
xlabel('\xi'), ylabel('y');
title("Non-dimensional temperature difference near tip by mesh size");
hold off;

% saves figure
saveas(gcf, 'question1.png');
