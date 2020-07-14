function [u_old, u_new, count] = successiveOverRelaxation(initial_guess, omega, epsilon)
    % initialize output variables
    u_old = initial_guess;
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
        
        % for all rows starting from bottom, except boundary
        for i = M-1:-1:2
            % for all columns starting at left, except boundary
            for j = 2:N-1
                % calculate new point
                u_new(i,j) = (u_old(i+1,j) + u_new(i-1,j) + u_old(i,j+1) + u_new(i,j-1))/4;
                u_new(i,j) = omega*u_new(i,j) + (1-omega)*u_old(i,j);
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