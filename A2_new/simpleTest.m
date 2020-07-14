% set constants
Bi = 0.25;
u_inf = 0;
epsilon = 1e-4;
length = 2;

initial_mesh = zeros(26, 51);

current = 200;
rate = 200/25;

initial_mesh(1,1) = 200;

for j = 1:26
    for i = 2:51
        initial_mesh(j,i) = current;
    end
    current = current - rate;
end

[u_old, u_new, count] = improvedBoundary(initial_mesh, Bi, u_inf, epsilon, length);
[M,N] = size(u_new);
figure;hold on;plot(linspace(0,1,M), u_new(:,M));plot(linspace(0,1,M), u_new(:,M),'*');