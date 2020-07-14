function [u_old, u_new, count] = alternatingDirectionImplicit(initial_guess, epsilon)
    % initialize output variables
    u_old = initial_guess;
    u_temp = initial_guess;
    u_new = initial_guess;
    count = 0;
    
    % Get mesh size from guess matrix
    [M, N] = size(initial_guess);
    
    % error variable declaration
    error = 100;
    
    while error>epsilon
        % set previously new solution as previous solution
        u_old = u_new;
        
        % increment iteration count
        count = count + 1;
        
        % Straight Line Relaxation on all Rows
        for j = 2:M-1   % for every row
            F = zeros(N,1);
            delta = zeros(N,1);
            for i = 2:N-1   % for every column, find F and delta
                F(i)= 1/(4-F(i-1));
                delta(i) = (u_old(j+1,i)+u_old(j-1,i)+delta(i-1))/(4-F(i-1));
            end
            
            for i = N-1:-1:2   % for every column, solve u
                u_temp(j, i) = u_temp(j, i+1)*F(i)+delta(i);
            end
        end
        
        % Straight Line Relaxation on all Columns
        for i = 2:N-1   % for every column
            F = zeros(M,1);
            delta = zeros(M,1);
            delta(1) = 200;
            for j = 2:M-1   % for every row, find F and delta
                F(j)= 1/(4-F(j-1));
                delta(j) = (u_temp(j,i+1)+u_temp(j,i-1)+delta(j-1))/(4-F(j-1));
            end
            
            for j = M-1:-1:2   % for every column, solve u
                u_new(j,i) = u_new(j+1,i)*F(j)+delta(j);
            end
        end
        
        % obtain error for element at (i=2, j=2)
        error = abs(1-u_old(2,2)/u_new(2,2));
        
        % calculate maximum difference between old and new solutions
        for i = 2:M-1
            for j = 2:N-1
                e_temp = abs(1-u_old(i,j)/u_new(i,j));
                if e_temp > error
                    error = e_temp;
                end
            end
        end
    end
    
    return;
end