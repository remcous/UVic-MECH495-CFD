%CONSTANTS and initial values
M = 200; % number of vertical & horizontal intervals M=N
psi_initial = ones(M+1,M+1);
zeta_initial = ones(M+1,M+1);
zeta_initial(:,1)=0;zeta_initial(:,M+1)=0;zeta_initial(1,:)=0;zeta_initial(M+1,:)=0;
k = 1/M;    % h=k
epsilon = 1e-4;
Re_start = 1;
Re_final = 13;

% boundary for psi = streamfunction
psi_initial(:,1) = 0; % set column 1 = 0
psi_initial(:,M+1) = 0; % set column M+1 = 0
psi_initial(1,:) = 0; % set row 1 = 0
psi_initial(M+1,:) = 0; % set row M+1 = 0

%initial boundary for zeta = vorticity transport
for j=2:M
    psi_initial(j,2) = psi_initial(j,3)/4;
    psi_initial(j,M) = psi_initial(j,M-1)/4;
    zeta_initial(j,2) = (4*psi_initial(j,2) - psi_initial(j,3) - psi_initial(j,1)...
        - psi_initial(j+1,2) - psi_initial(j-1,2)) / (k.^2);
    zeta_initial(j,M) = (4*psi_initial(j,M) - psi_initial(j,M+1) - psi_initial(j,M-1)...
        - psi_initial(j+1,M) - psi_initial(j-1,M)) / (k.^2);
end
for i=2:M
    psi_initial(2,i) = psi_initial(3,i)/4;
    psi_initial(M,i) = psi_initial(M-1,i)/4 - k/2;
    zeta_initial(2,i) = (4*psi_initial(2,i) - psi_initial(3,i) - psi_initial(1,i)...
        - psi_initial(2,i+1) - psi_initial(2,i-1)) / (k.^2);
    zeta_initial(M,i) = (4*psi_initial(M,i) - psi_initial(M+1,i) - psi_initial(M-1,i)...
        - psi_initial(M,i+1) - psi_initial(M,i-1)) / (k.^2);
end

drivenInterior(psi_initial, zeta_initial, Re_start, Re_final, epsilon);