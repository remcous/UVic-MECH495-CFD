function [u_new, count] = explicit(u_0, delta_t, t_final)
    % get necessary information for calculations
    u_new = u_0;
    M = size(u_0);
    M = M(1) - 1;
    delta_x = 1/M;
    
    % initialize error and iteration count
    count = 0;
    plot(linspace(0,1,M+1),u_0); hold on;
    
    for t = 0:delta_t:t_final
        % set previous iteration matrix and increment iteration count
        u_old = u_new;
        count = count + 1;
        
        % calculate solution for new time step at each point
        for i = 2:M
            u_new(i) = delta_t/(delta_x.^2) * (u_old(i+1) - 2*u_old(i) + u_old(i-1)) + u_old(i);
        end
        
        if mod(t,0.1) == 0
            plot(linspace(0,1,M+1),u_new);
        end
    end
    
    title(['First Order explicit solution with \Deltat = ' num2str(delta_t)]);
    legend('t=0', 't=0.1', 't=0.2', 't=0.3', 't=0.4', 't=0.5', 't=0.6', ...
        't=0.7', 't=0.8', 't=0.9', 't=1.0');
    xlabel('x'); ylabel('u(x,t)'); hold off;
    saveas(gcf,[pwd '\explicit_' num2str(delta_t*1000) '.png']); close(gcf);
end