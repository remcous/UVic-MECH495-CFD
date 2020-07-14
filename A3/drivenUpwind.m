function drivenUpwind(psi_0,zeta_0,Re_i, Re_f, epsilon)
    % get size square cavity
    [N,M] = size(psi_0);
    M = M -1;   % N=M
    k = 1/M;    % h=k
    u_new = zeros(M+1,M+1); u_new(M+1,:) = 1;
    v_new = zeros(M+1,M+1);
    C = zeros(Re_f-Re_i,1);
    
    for Re = Re_i:Re_f
        count =0;
        error = 100;
        psi_new = psi_0;
        zeta_new = zeta_0;
        
        while(error>epsilon)
            % set previous iteration to old
            psi_old = psi_new;
            zeta_old = zeta_new;
            u_old = u_new;
            v_old = v_new;
            
            if count == 5e5
                fprintf(['The iterations have exceeded 5e5 for a Re of ' num2str(Re) '\nBreaking out of loop']);
                break;
            end

            % increment iteration counter
            count = count + 1;

            %calculate boundary for zeta = vorticity transport
            for j=1:M+1
                zeta_new(j,1) = -2*psi_old(j,2)/(k.^2);
                zeta_new(j,M+1) = -2*psi_old(j,M)/(k.^2);
            end
            for i=1:M+1
                zeta_new(1,i) = -2*psi_old(2,i)/(k.^2);
                zeta_new(M+1,i) = -2*(k+psi_old(M,i)+psi_old(M+1,i))/(k.^2);
            end

            % calculate values of zeta
            for j=2:M % rows
                for i=2:M % columns
                    u_term = (psi_old(j+1,i)-psi_old(j-1,i))/2;
                    v_term = (psi_old(j,i-1)-psi_old(j,i+1))/2;
                    
                    if u_term > 0
                        if v_term > 0
                            zeta_new(j,i) = ...
                                (u_term * zeta_new(j,i-1) + v_term * zeta_new(j-1,i)...
                                + 1/Re * (zeta_old(j,i+1)+zeta_new(j,i-1)+...
                                zeta_old(j+1,i)+zeta_new(j-1,i))) / (u_term + v_term + 4/Re);
                        else
                            zeta_new(j,i) = ...
                                (u_term * zeta_new(j,i-1) - v_term * zeta_new(j+1,i)...
                                + 1/Re * (zeta_old(j,i+1)+zeta_new(j,i-1)+...
                                zeta_old(j+1,i)+zeta_new(j-1,i))) / (u_term - v_term + 4/Re);
                        end
                    else
                        if v_term > 0
                            zeta_new(j,i) = ...
                                (- u_term * zeta_new(j,i+1) + v_term * zeta_new(j-1,i)...
                                + 1/Re * (zeta_old(j,i+1)+zeta_new(j,i-1)+...
                                zeta_old(j+1,i)+zeta_new(j-1,i))) / (- u_term + v_term + 4/Re);
                        else
                            zeta_new(j,i) = ...
                                (- u_term * zeta_new(j,i+1) - v_term * zeta_new(j+1,i)...
                                + 1/Re * (zeta_old(j,i+1)+zeta_new(j,i-1)+...
                                zeta_old(j+1,i)+zeta_new(j-1,i))) / (- u_term - v_term + 4/Re);
                        end
                    end
                end
            end

            % calculate values of psi
            for j=2:M % rows
                for i=2:M % columns
                    psi_new(j,i) = ...
                        (psi_old(j,i+1)+psi_old(j,i-1)+psi_old(j+1,i)+psi_old(j-1,i)+(k.^2)*zeta_new(j,i)) / 4;
                end
            end

            % initialize error at 0
            error = 0;

            % Calculate u,v for error calculation
            for j= 2:M
                for i = 2:M
                    u_new(j,i)=(psi_new(j+1,i)-psi_new(j-1,i))/(2*k);
                    v_new(j,i)=-(psi_new(j,i+1)-psi_new(j,i-1))/(2*k);
                end
            end
        
            % calculate maximum difference between old and new solutions
            for i = 1:M+1
                for j = 1:M+1
                    % check u
                    e_temp = abs(1-u_old(j,i)/u_new(j,i));
                    if e_temp > error
                        error = e_temp;
                    end
                    % check v
                    e_temp = abs(1-v_old(j,i)/v_new(j,i));
                    if e_temp > error
                        error = e_temp;
                    end
                end
            end
        end
        
        C(Re-Re_i+1,1) = count;
        
        % Create and save figure for current Reynold's number
        figure;contourf(linspace(0,1,M+1), linspace(0,1,M+1), zeta_new);colorbar;
        daspect([1 1 1]);set(gcf, 'Position', get(0, 'Screensize'));
        title(['\zeta for Re=', num2str(Re), ' after ', num2str(count), ' iterations using Upwind Method']);
        saveas(gcf,[pwd '\Figures\Upwind\zeta\zeta_' num2str(Re) '.png']); close(gcf);
        figure;contourf(linspace(0,1,M+1), linspace(0,1,M+1), psi_new);colorbar;
        daspect([1 1 1]);set(gcf, 'Position', get(0, 'Screensize'));
        title(['\Psi for Re=', num2str(Re), ' after ', num2str(count), ' iterations using Upwind Method']);
        saveas(gcf,[pwd '\Figures\Upwind\psi\psi_' num2str(Re) '.png']); close(gcf);
        figure;q=quiver(linspace(0,1,M+1),linspace(0,1,M+1),u_new,v_new);colorbar;hold on;
        q.AutoScaleFactor = 4;
        contour(linspace(0,1,M+1), linspace(0,1,M+1), psi_new);hold off;
        daspect([1 1 1]);set(gcf, 'Position', get(0, 'Screensize'));
        xlim([0 1]); ylim([0 1]);
        title(['Vector Field for Re=', num2str(Re), ' after ', num2str(count), ' iterations using Upwind Method']);
        saveas(gcf,[pwd '\Figures\Upwind\vector\vector_' num2str(Re) '.png']); close(gcf);
        
    end
    
    if Re_i ~= Re_f
        plot(linspace(Re_i,Re_f,Re_f-Re_i+1),C);
        xlim([Re_i Re_f]);xlabel('Re'); ylabel('# of Iterations');
        title('Iterations to Converge by Reynolds Number');
    end
end