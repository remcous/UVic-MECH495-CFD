% acceptable error for solution
epsilon = 1e-4;

% chosen value of omega for SOR
omega = 1.32;

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

% 50x25 intervals plots Jakobi
[u_old, u_new, count] = jakobi(initial_mesh_25, epsilon);
[M,N]=size(u_new);
title_prefix = 'Jacobi '; folder='Jakobi'; prefix='jakobi'; suffix=25;
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
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ...
    ' mesh intervals difference for iteration ', num2str(count-1), ...
    ' to ' num2str(count)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_diff.png']); close(gcf);

% 100x50 intervals plots Jakobi
[u_old, u_new, count] = jakobi(initial_mesh_50, epsilon);
[M,N]=size(u_new);
title_prefix = 'Jacobi '; folder='Jakobi'; prefix='jakobi'; suffix=50;
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
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ...
    ' mesh intervals difference for iteration ', num2str(count-1), ...
    ' to ' num2str(count)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_diff.png']); close(gcf);

% 200x100 intervals plots Jakobi
[u_old, u_new, count] = jakobi(initial_mesh_100, epsilon);
[M,N]=size(u_new);
title_prefix = 'Jacobi '; folder='Jakobi'; prefix='jakobi'; suffix=100;
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
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ...
    ' mesh intervals difference for iteration ', num2str(count-1), ...
    ' to ' num2str(count)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_diff.png']); close(gcf);

% 50x25 intervals plots Gauss-Seidel
[u_old, u_new, count] = gaussSeidel(initial_mesh_25, epsilon);
[M,N]=size(u_new);
title_prefix = 'Gauss-Seidel '; folder='GaussSeidel'; prefix='gs'; suffix=25;
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
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ...
    ' mesh intervals difference for iteration ', num2str(count-1), ...
    ' to ' num2str(count)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_diff.png']); close(gcf);

% 100x50 intervals plots Gauss-Seidel
[u_old, u_new, count] = gaussSeidel(initial_mesh_50, epsilon);
[M,N]=size(u_new);
title_prefix = 'Gauss-Seidel '; folder='GaussSeidel'; prefix='gs'; suffix=50;
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
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ...
    ' mesh intervals difference for iteration ', num2str(count-1), ...
    ' to ' num2str(count)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_diff.png']); close(gcf);

% 200x100 intervals plots Gauss-Seidel
[u_old, u_new, count] = gaussSeidel(initial_mesh_100, epsilon);
[M,N]=size(u_new);
title_prefix = 'Gauss-Seidel '; folder='GaussSeidel'; prefix='gs'; suffix=100;
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
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ...
    ' mesh intervals difference for iteration ', num2str(count-1), ...
    ' to ' num2str(count)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_diff.png']); close(gcf);


% 50x25 intervals plots Successive Over Relaxation
[u_old, u_new, count] = successiveOverRelaxation(initial_mesh_25, omega, epsilon);
[M,N]=size(u_new);
title_prefix = 'Successive Over Relaxation '; folder='SOR'; prefix='sor'; suffix=25;
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
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ...
    ' mesh intervals difference for iteration ', num2str(count-1), ...
    ' to ' num2str(count)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_diff.png']); close(gcf);

% 100x50 intervals plots Successive Over Relaxation
[u_old, u_new, count] = successiveOverRelaxation(initial_mesh_50, omega, epsilon);
[M,N]=size(u_new);
title_prefix = 'Successive Over Relaxation '; folder='SOR'; prefix='sor'; suffix=50;
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
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ...
    ' mesh intervals difference for iteration ', num2str(count-1), ...
    ' to ' num2str(count)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_diff.png']); close(gcf);

% 200x100 intervals plots Successive Over Relaxation
[u_old, u_new, count] = successiveOverRelaxation(initial_mesh_100, omega, epsilon);
[M,N]=size(u_new);
title_prefix = 'Successive Over Relaxation '; folder='SOR'; prefix='sor'; suffix=100;
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
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ...
    ' mesh intervals difference for iteration ', num2str(count-1), ...
    ' to ' num2str(count)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_diff.png']); close(gcf);


% 50x25 intervals plots Alternating Direction Implicit Method
[u_old, u_new, count] = alternatingDirectionImplicit(initial_mesh_25, epsilon);
[M,N]=size(u_new);
title_prefix = 'Alternating Direction Implicit '; folder='ADI'; prefix='adi'; suffix=25;
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
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ...
    ' mesh intervals difference for iteration ', num2str(count-1), ...
    ' to ' num2str(count)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_diff.png']); close(gcf);

% 100x50 intervals plots Alternating Direction Implicit Method
[u_old, u_new, count] = alternatingDirectionImplicit(initial_mesh_50, epsilon);
[M,N]=size(u_new);
title_prefix = 'Alternating Direction Implicit '; folder='ADI'; prefix='adi'; suffix=50;
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
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ...
    ' mesh intervals difference for iteration ', num2str(count-1), ...
    ' to ' num2str(count)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_diff.png']); close(gcf);

% 200x100 intervals plots Alternating Direction Implicit Method
[u_old, u_new, count] = alternatingDirectionImplicit(initial_mesh_100, epsilon);
[M,N]=size(u_new);
title_prefix = 'Alternating Direction Implicit '; folder='ADI'; prefix='adi'; suffix=100;
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
title([title_prefix, num2str(N-1), 'x', num2str(M-1), ...
    ' mesh intervals difference for iteration ', num2str(count-1), ...
    ' to ' num2str(count)]);
saveas(gcf,[pwd '\' folder '\' prefix num2str(suffix) '_diff.png']); close(gcf);