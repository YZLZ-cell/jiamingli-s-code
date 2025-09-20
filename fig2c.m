%%Parameter setup
clear;
clc;
w = linspace(1, 3, 20);
d = linspace(-20, 20, 40);
t1 = 1;
t2 = 4;

ave = zeros(length(w), length(d));
y_val1_arr = zeros(1, length(w));
y_val2_arr = zeros(1, length(w));

for i = 1:length(w)
    w_val = w(i);
    N = round(1000 * (t2 + 2*t1));
    T = (0.4 - (w_val-2)*0.05);
    r = 1.1*w_val*w_val/4;
    r1 = w_val*w_val/4*0.6;
    exp_term = exp(-t1/T);
    exp_t1_over_T = exp(t1/T);
    d0 = r1 * t2/T - r*2*(exp_term - exp_t1_over_T)/(exp_term + exp_t1_over_T);
    for j = 1:length(d)
        d_val = d(j);
        ave(i, j) = eve(N, T, t1, t2, r, r1, w_val, d_val, 0);
    end
end

figure;
hold on;

contourf(w, d, ave', ...
    'LevelList', linspace(0, 1, 11), 'LineColor', 'none');

ax = gca;
ax.XLabel.String = '$\Omega$ (MHz)';
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = 25;
ax.XLabel.Position(2) = ax.XLabel.Position(2)+3.5 ;

ax.YLabel.String = '$\delta$ (MHz)';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = 25;
ax.YLabel.Position(1) = ax.YLabel.Position(1)+0.1; 

cb = colorbar;
caxis([0 1]);
cb.Ticks = 0:0.5:1;
cb.Label.FontSize = 25;

ax.XTick = linspace(1, 3, 2);
ax.YTick = linspace(-20, 20, 3);
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 25;

drawnow;

%%
function out = sig_0()
    out = diag([1, 1]);
end

function out = sig_z()
    out = diag([1, -1])/2;
end

function out = sig_up()
    out = [0, 1; 0, 0]/2; 
end

function out = sig_down()
    out = [0, 0; 1, 0]/2;
end

function out = com(a, b)
    out = a*b - b*a;
end

function h = ham(t, T, t1, t2, r, r1, w, d, d0)
    a = r * T * log(1 / (exp((t - t1) / T) + exp(-(t - t1) / T))) + r1 / T * t2 / 2 * t - ...
        r * T * log(1 / (exp((0 - t1) / T) + exp(-(0 - t1) / T)));
    b = r * T * log(1 / (exp((t - t1 - t2) / T) + exp(-(t - t1 - t2) / T))) - ...
        r1 / T * t2 / 2 * (t - t2 - t1);
    c = -r1 / T * (t - t1 - t2 / 2)^2 / 2 + r1 / T * (-t2 / 2)^2 / 2;
    e = r * T * log(1 / (exp((t1 - t1) / T) + exp(-(t1 - t1) / T))) + ...
        r1 / T * t2 / 2 * t1 - r * T * log(1 / (exp((0 - t1) / T) + exp(-(0 - t1) / T)));

    if t < t1
        x = 2 * w / (exp((t - t1) / T) + exp(-(t - t1) / T)) * (exp(-1i * a + 1i * d0 / 2 * t)) ;
        y = 2 * w / (exp((t - t1) / T) + exp(-(t - t1) / T)) * (exp(1i * a - 1i * d0 / 2 * t)) ;
        z = d;
    elseif t1 <= t && t < t1 + t2
        x = w * (exp(-1i * (c + e) + 1i * d0 / 2 * t)) ;
        y = w * (exp(1i * (c + e) - 1i * d0 / 2 * t)) ;
        z = d;
    else
        x = 2 * w / (exp((t - t1 - t2) / T) + exp(-(t - t1 - t2) / T)) * (exp(-1i * (b + e) + 1i * d0 / 2 * t));
        y = 2 * w / (exp((t - t1 - t2) / T) + exp(-(t - t1 - t2) / T)) * (exp(1i * (b + e) - 1i * d0 / 2 * t) );
        z = d;
    end

    h = x * sig_up() + y * sig_down() + z * sig_z();
end

function result = eve(N, T, t1, t2, r, r1, w, d, d0)
    rho0 = sig_z() + sig_0() / 2;
    rho = rho0;
    result_arr = zeros(N, 1);

    for i = 1:N
        dt = 0.001;
        h = ham(i * dt, T, t1, t2, r, r1, w, d, d0);
        rho = rho - 1i * com(h, rho) * dt - 0.5 * com(h, com(h, rho)) * (dt^2);
        result_arr(i) = real(rho(2, 2));
    end

    result = result_arr(end);
end