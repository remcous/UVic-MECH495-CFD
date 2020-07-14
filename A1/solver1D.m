function [xi, y] = solver1D(a, b, c, Bi, N)
    xi = linspace(0,1,N+1);
    d_xi = xi(2);

    % initialize F and delta vectors
    F = zeros(N,1);
    delta = zeros(N,1);
    delta(1) = 1;

    % Calculate all F and delta values
    for j = 2:N
        if j == N
            F(j) = -b(xi(j-1),d_xi)/(a(xi(j-1),d_xi)+c(xi(j-1),d_xi)*F(j-1));
            delta(j) = -c(xi(j-1),d_xi)*delta(j-1)/(a(xi(j-1),d_xi)+c(xi(j-1),d_xi)*F(j-1));
        else
            F(j) = -b(xi(j),d_xi)/(a(xi(j),d_xi)+c(xi(j),d_xi)*F(j-1));
            delta(j) = -c(xi(j),d_xi)*delta(j-1)/(a(xi(j),d_xi)+c(xi(j),d_xi)*F(j-1));
        end
    end

    % establish y vector and boundary conditions
    y = zeros(N+1,1);
    y(1)=1;
    num = (3+F(N-1)*(1/3*F(N)-3/2))*delta(N)+(1/3*F(N-2)-3/2)*delta(N-1)+1/3*delta(N-2);
    denom = 11/6+Bi*d_xi+F(N)*(-3+F(N-1)*(3/2-1/3*F(N-2)));
    y(N+1) = num/denom;

    % calculate y values
    for j = N:-1:2
        y(j) = y(j+1)*F(j)+delta(j);
    end
end