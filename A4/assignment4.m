% constants
delta_t = 0.005;
delta_x = 0.1;
t_final = 1;

u_0 = zeros(1/delta_x+1,1);
for i=1:size(u_0)
    x = delta_x * (i-1);
    if x <= 0.5
        u_0(i) = 2*x;
    else
        u_0(i) = 2*(1-x);
    end
end

explicit(u_0, delta_t, t_final);
implicit(u_0, delta_t, t_final);
explicit(u_0, 0.01, t_final);
implicit(u_0, 0.01, t_final);