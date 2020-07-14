%CONSTANTS and initial values
M = 100; % number of vertical & horizontal intervals M=N
k = 1/M;    % h=k
epsilon = 1e-4;
Re_start = 1;
Re_final = 15;
psi_initial = ones(M+1,M+1);
zeta_initial = ones(M+1,M+1);

% boundary for psi = streamfunction
psi_initial(:,1) = 0; % set column 1 = 0
psi_initial(:,M+1) = 0; % set column M+1 = 0
psi_initial(1,:) = 0; % set row 1 = 0
psi_initial(M+1,:) = 0; % set row M+1 = 0

%initial boundary for zeta = vorticity transport
for j=1:M+1
    zeta_initial(j,1) = -2*psi_initial(j,2)/(k.^2);
    zeta_initial(j,M+1) = -2*psi_initial(j,M)/(k.^2);
end
for i=1:M+1
    zeta_initial(1,i) = -2*psi_initial(2,i)/(k.^2);
    zeta_initial(M+1,i) = -2*(k+psi_initial(M,i)+psi_initial(M+1,i))/(k.^2);
end

drivenCavity(psi_initial, zeta_initial, Re_start, Re_final, epsilon);