function [u] = implicit(u_0, delta_t, t_final)
    % get necessary information for calculations
    u = u_0;
    M = size(u_0);
    M = M(1) - 1;
    delta_x = 1/M;
    
    % initialize F and delta arrays (F(1) and delta(1) = 0)
    F = zeros(M,1);
    delta = F;
    a = 1 + delta_t/(delta_x.^2);
    b = -delta_t/(2 * delta_x.^2);
    c = b;
    d = zeros(M+1,1); % a,b,c don't change with i, however d does
    
    % Plot the initial conditions
    plot(linspace(0,1,M+1),u_0); hold on;
    
    for t = 0:delta_t:t_final
        % calculate the values required for Thomas algorithm
        for i = 2:M
            d(i) = delta_t / (2 * delta_x.^2) * (u(i+1) - 2*u(i)...
                + u(i-1)) + u(i);
            F(i) = delta_t / (2*delta_x.^2 + 2*delta_t*(1-0.5 * F(i-1)));
            delta(i) = (delta_t*(u(i+1)-2*u(i)+u(i-1))...
                + 2*delta_x.^2*u(i)+delta_t*delta(i-1))...
                / (2*delta_x.^2 + 2*delta_t*(1-0.5 * F(i-1)));
        end
        
        % Calculate new values of u for all i
        for i = M:-1:2
            u(i) = u(i+1)*F(i)+delta(i);
        end
        
        % plot multiples of 0.1 in time
        if mod(t,0.1) == 0
            plot(linspace(0,1,M+1),u);
        end
    end
    
    title(['Crank-Nicholson Implicit solution with \Deltat = ' num2str(delta_t)]);
    legend('t=0', 't=0.1', 't=0.2', 't=0.3', 't=0.4', 't=0.5', 't=0.6', ...
        't=0.7', 't=0.8', 't=0.9', 't=1.0');
    xlabel('x'); ylabel('u(x,t)'); hold off;
    saveas(gcf,[pwd '\implicit_' num2str(delta_t*1000) '.png']); close(gcf);
end