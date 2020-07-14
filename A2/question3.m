% set constants
Bi = 0.25;
u_inf = 0;
epsilon = 1e-8;
length = 2;

% produce matrices for initial guesses for each interval size
initial_mesh_100 = zeros(101,201);
initial_mesh_50 = zeros(51,101);
initial_mesh_25 = zeros(26, 51);

% populate the matrices with the boundary conditions
for j = 1:51
    initial_mesh_25(1,j) = 200;
    initial_mesh_50(1,j) = 200;
    initial_mesh_100(1,j) = 200;
end
for j = 52:101
    initial_mesh_50(1,j) = 200;
    initial_mesh_100(1,j) = 200;
end
for j = 102:201
    initial_mesh_100(1,j) = 200;
end

% set central elements to a constant
for i = 2:25
    for j = 2:50
        initial_mesh_25(i,j) = 100;
    end
end
for i = 2:50
    for j = 2:100
        initial_mesh_50(i,j) = 100;
    end
end
for i = 2:100
    for j = 2:200
        initial_mesh_100(i,j) = 100;
    end
end

% 50x25 intervals plots Improved Boundary Treatment
[u_old, u_new, count] = improvedBoundary(initial_mesh_25, Bi, u_inf, epsilon, length);
[M,N]=size(u_new);
title_prefix = 'Improved Boundary Treatment '; folder='Improved'; 
prefix='improved'; suffix=25;
figure;pcolor(linspace(0,2,N), linspace(0,1,M), u_new);colorbar;
daspect([1 1 1]);set(gcf, 'Position', get(0, 'Screensize'));
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ' mesh intervals at iteration ', num2str(count)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_new.png']); close(gcf);
figure;pcolor(linspace(0,2,N), linspace(0,1,M), u_old);colorbar;
daspect([1 1 1]);set(gcf, 'Position', get(0, 'Screensize'));
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ' mesh intervals at iteration ', num2str(count-1)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_old.png']); close(gcf);
figure;pcolor(linspace(0,2,N), linspace(0,1,M), u_new-u_old);colorbar;
daspect([1 1 1]);set(gcf, 'Position', get(0, 'Screensize'));
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ' mesh intervals final difference']);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_diff.png']); close(gcf);

% 100x50 intervals plots Improved Boundary Treatment
[u_old, u_new, count] = improvedBoundary(initial_mesh_50, Bi, u_inf, epsilon, length);
[M,N]=size(u_new);
title_prefix = 'Improved Boundary Treatment '; folder='Improved'; 
prefix='improved'; suffix=50;
figure;pcolor(linspace(0,2,N), linspace(0,1,M), u_new);colorbar;
daspect([1 1 1]);set(gcf, 'Position', get(0, 'Screensize'));
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ' mesh intervals at iteration ', num2str(count)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_new.png']); close(gcf);
figure;pcolor(linspace(0,2,N), linspace(0,1,M), u_old);colorbar;
daspect([1 1 1]);set(gcf, 'Position', get(0, 'Screensize'));
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ' mesh intervals at iteration ', num2str(count-1)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_old.png']); close(gcf);
figure;pcolor(linspace(0,2,N), linspace(0,1,M), u_new-u_old);colorbar;
daspect([1 1 1]);set(gcf, 'Position', get(0, 'Screensize'));
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ' mesh intervals final difference']);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_diff.png']); close(gcf);

% 200x100 intervals plots Improved Boundary Treatment
[u_old, u_new, count] = improvedBoundary(initial_mesh_100, Bi, u_inf, epsilon, length);
[M,N]=size(u_new);
title_prefix = 'Improved Boundary Treatment '; folder='Improved'; 
prefix='improved'; suffix=100;
figure;pcolor(linspace(0,2,N), linspace(0,1,M), u_new);colorbar;
daspect([1 1 1]);set(gcf, 'Position', get(0, 'Screensize'));
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ' mesh intervals at iteration ', num2str(count)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_new.png']); close(gcf);
figure;pcolor(linspace(0,2,N), linspace(0,1,M), u_old);colorbar;
daspect([1 1 1]);set(gcf, 'Position', get(0, 'Screensize'));
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ' mesh intervals at iteration ', num2str(count-1)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_old.png']); close(gcf);
figure;pcolor(linspace(0,2,N), linspace(0,1,M), u_new-u_old);colorbar;
daspect([1 1 1]);set(gcf, 'Position', get(0, 'Screensize'));
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ' mesh intervals final difference']);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_diff.png']); close(gcf);