%Parameter setup
clear;
clc;

t1 = 0.5;
t2 = linspace(2, 6, 8);
B = 80;
linspace(0.01, 0.5, 6);
w = linspace(2, 5, 8);
p = linspace(pi, 2*pi, 10);
d = linspace(0, B/2, 10);



ave = zeros(length(w), length(t2));
av = zeros(length(w), length(t2), 4);
avd = zeros(length(w), length(t2), 4, 5);
current_avg = zeros(length(w), length(t2), 4, 5, length(p));
final_transfer = zeros(length(w), length(t2), 4,  5, length(p), length(d));

for i = 1:length(w)
    w_val = w(i);
    r = linspace(w_val-1, w_val+0.5, 5);
    T = linspace(0.2-(w_val-2)*0.05, 0.4-(w_val-2)*0.05, 4 );
    for j = 1:length(t2)
        t2_val = t2(j);
        N = round(1000 * (t2_val + 2 * t1));
        for n = 1:length(T)
            T_val = T(n);
            for m = 1:length(r)
                r_val = r(m);
                for l = 1:length(p)
                    p_val = p(l);
                    for k = 1:length(d)
                        d_val = d(k); 
                        N = round(1000 * (t2_val + 2 * t1));
                        d0_val = B/3.2 - r_val * 2 * (exp(-t1/T_val) - exp(-(-t1)/T_val)) / (exp(-t1/T_val) + exp(-(-t1)/T_val));
                        final_transfer(i, j, n, m, l, k) = eve(N, T_val, t1, t2_val,  r_val , B/3.2 * T_val / t2_val, w_val, d_val, d0_val, p_val);
                    end
                    current_avg(i, j, n, m, l) = mean(final_transfer(i, j, n, m, l, :));
                end
                avd(i, j, n, m) = mean(current_avg(i, j, n, m, :));
            end
            av(i, j, n) = max(avd(i, j, n, :));
        end
        ave(i, j) = max(av(i, j, :));
    end
end


figure;
contourf( w, t2+2*t1, ave.', 10, 'LineColor', 'none');
ax = gca;
box off;

ax.XLabel.String = '$\Omega$ (MHz)';
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = 25;
ax.XLabel.Position(2) = ax.XLabel.Position(2)+0.35 ;

ax.YLabel.String = '$\tau (\mu s)$';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize =25;
ax.YLabel.Position(1) = ax.YLabel.Position(1)+0.05; 


cb = colorbar;
caxis([0 1]);
cb.Ticks = 0:0.5:1;
cb.Label.FontSize = 25;

ax.YTick = linspace(3, 7, 2);
ax.XTick = linspace(2, 5, 2);
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 25;

drawnow;
%%
function out = sig_0()
    out = diag([1, 1]);
end

function out = sig_z()
    out = diag([1, -1]) / 2;
end

function out = sig_x()
    out = [0, 1; 1, 0] / 2;
end

function out = sig_up()
    out = [0, 1; 0, 0] / 2;
end

function out = sig_down()
    out = [0, 0; 1, 0] / 2;
end

function out = com(a, b)
    out = a * b - b * a;
end

function h3 = ham3(t, T, t1, t2, r, r1, w, d, d0, p)
    a = r*T*log(1/(exp((t-t1)/T) + exp(-(t-t1)/T))) + r1/T*t2/2*t - r*T*log(1/(exp((0-t1)/T) + exp(-(0-t1)/T)));
    b = r*T*log(1/(exp((t-t1-t2)/T) + exp(-(t-t1-t2)/T))) - r1/T*t2/2*(t-t2-t1);
    c = -r1/T*(t-t1-t2/2)^2/2 + r1/T*(-t2/2)^2/2;
    e = r*T*log(1/(exp((t1-t1)/T) + exp(-(t1-t1)/T))) + r1/T*t2/2*t1 - r*T*log(1/(exp((0-t1)/T) + exp(-(0-t1)/T)));
    
    if t < t1
        x = 2*w/(exp((t-t1)/T) + exp(-(t-t1)/T)) * (exp(1i*a+1i*p  ) + 0.95*exp(-1i*a - 1i*d0*t) + 0.9*exp(-1i*a + 1i*(d0)*t) );
        y = 2*w/(exp((t-t1)/T) + exp(-(t-t1)/T)) * (exp(-1i*a-1i*p ) + 0.95*exp(1i*a + 1i*d0*t) + 0.9*exp(1i*a - 1i*(d0)*t) );
        z = d;
    elseif t1 <= t && t < t1 + t2
        x = w * (exp(1i*(c + e)+1i*p) + 0.95*exp(-1i*(c + e) - 1i*d0*t) + 0.9*exp(-1i*(c + e) + 1i*d0*t) );
        y = w * (exp(-1i*(c + e)-1i*p ) + 0.95*exp(1i*(c + e) + 1i*d0*t) + 0.9*exp(1i*(c + e) - 1i*d0*t) );
        z = d;
    else
        x = 2*w/(exp((t-t1-t2)/T) + exp(-(t-t1-t2)/T)) * (exp(1i*(b + e)+1i*p) + 0.95*exp(-1i*(b + e) - 1i*d0*t) + 0.9*exp(-1i*(b + e) + 1i*d0*t) );
        y = 2*w/(exp((t-t1-t2)/T) + exp(-(t-t1-t2)/T)) * (exp(-1i*(b + e)-1i*p ) + 0.95*exp(1i*(b + e) + 1i*d0*t) + 0.9*exp(1i*(b + e) - 1i*d0*t) );
        z = d;
    end
    
    h3 = x*sig_up() + y*sig_down() + z*sig_z() - 1i*0.1*(sig_z()+sig_0())/2;
end
function result3 = eve(N, T, t1, t2, r, r1, w, d, d0, p)
    rho0 = sig_z() + sig_0() / 2;
    rho = rho0;
    result_arr3 = zeros(N, 1);

    for i = 1:N
        dt = 0.001;
        h3 = ham3(i * dt, T, t1, t2, r, r1, w, d, d0, p);
        rho = rho - 1i * com(h3, rho) * dt - 0.5 * com(h3, com(h3, rho)) * (dt^2);
        result_arr3(i) = real(rho(2, 2));
    end

    result3 = result_arr3(end);
end