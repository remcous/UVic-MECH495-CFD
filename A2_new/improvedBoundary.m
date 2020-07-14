function [u_old, u_new, count] = improvedBoundary(u_initial, Bi, ...
    u_inf, epsilon, length)
    % get mesh size from initial guess
    [M,N] = size(u_initial); % N rows, M columns
    h = length/N;
    
    count = 0;
    
    % create matrices to hold previous and current solution
    u_old = u_initial;
    u_new = u_initial;
    
    % create error variable
    error = 100;
    
    while error > epsilon
        % store previous iteration as u_old
        u_old = u_new;
        
        % increment count
        count = count + 1;
        
        % Calculate boundary at j = N
        for i = 2:M-1
            u_new(i,N) = (2*u_old(i,N-3) - 9*u_old(i,N-2) ...
                + 18*u_old(i,N-1)) / (11+6*h*Bi);
        end
        
        % Calculate boundary at j = N-1
        for i = 2:M-1
            u_new(i,N-1) = (u_new(i,N) + u_old(i,N-2) + u_old(i+1,N-1) ...
                + u_new(i-1,N-1)) / 4;
        end
        
        % for all rows starting from bottom, except boundary
        for i = M-1:-1:2
            % for all columns starting at left, except boundary
            for j = 2:N-2
                % calculate new point
                u_new(i,j) = (u_old(i+1,j) + u_new(i-1,j) ...
                    + u_old(i,j+1) + u_new(i,j-1))/4;
            end
        end
        
        % obtain error for element at (i=2, j=2)
        error = abs(1-u_old(2,2)/u_new(2,2));
        
        % check error at all points
        for i = 2:M-1
            for j = 2:N-1
                e_temp = abs(1-u_old(i,j)/u_new(i,j));
                if e_temp > error
                    error = e_temp;
                end
            end
        end    
    end
    
    return
end