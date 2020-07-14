function [y_new, error] = solver1Dnonlinear(a, b, c, alpha, xi, y_old, N, epsilon)
    d_xi = xi(2)-xi(1);
    
    % initialize F and delta vectors
    F = zeros(N,1);
    delta = zeros(N,1);
    delta(1) = 1;
    
    error = 100;
    
    while(error>epsilon)
        % Calculate all F and delta values
        for j = 2:N
            if j == N
                F(j) = -b(xi(j-1),d_xi)/(a(xi(j-1),d_xi, alpha, y_old(j))+c(xi(j-1),d_xi)*F(j-1));
                delta(j) = -c(xi(j-1),d_xi)*delta(j-1)/(a(xi(j-1),d_xi, alpha, y_old(j))+c(xi(j-1),d_xi)*F(j-1));
            else
                F(j) = -b(xi(j),d_xi)/(a(xi(j),d_xi, alpha, y_old(j))+c(xi(j),d_xi)*F(j-1));
                delta(j) = -c(xi(j),d_xi)*delta(j-1)/(a(xi(j),d_xi, alpha, y_old(j))+c(xi(j),d_xi)*F(j-1));
            end
        end
        
        % establish y vector and boundary conditions
        y_new = zeros(N+1,1);
        y_new(1)=1;
        num = (3+F(N-1)*(1/3*F(N)-3/2))*delta(N)+(1/3*F(N-2)-3/2)*delta(N-1)+1/3*delta(N-2);
        denom = 11/6+(1/(1+alpha*y_old))*d_xi+F(N)*(-3+F(N-1)*(3/2-1/3*F(N-2)));
        %y_new(N+1) = num/denom;
    
        % calculate y values
        for j = N:-1:2
            y_new(j) = y_new(j+1)*F(j)+delta(j);
        end
        
        error = abs(1-y_old(2)/y_new(2));
        
        % find the minimum error 
        for j = 3:N
            e_temp = abs(1-y_old(j)/y_new(j));
            if e_temp > error
                error = e_temp;
            end
        end
        
        % update y_old
        y_old = y_new;
    end
end