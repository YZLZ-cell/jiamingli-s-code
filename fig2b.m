%Parameter setup
clear;
clc;
w = 2;
r = 1.2;
d_values = linspace(-20, 20, 40);
T = 0.35;
t1 = 1;
t2 = linspace(1, 7, 40);
r1 = 0.7;

average_population_transfers = zeros(length(t2), length(d_values));

for i = 1:length(t2)
    t2_val = t2(i);
    N = round(1000 * (t2_val + 2*t1));
    for j = 1:length(d_values)
        d_val = d_values(j);
        exp_term = exp(-t1/T);
        d0 = 0;
        average_population_transfers(i, j) = eve(N, T, t1, t2_val, r, r1, w, d_val, d0);
    end
end


figure;
contourf(t2+2*t1, d_values, average_population_transfers.', ...
    'LevelList', linspace(0, 1, 11), 'LineColor', 'none');

ax = gca;

box off;
ax.XLabel.String = '$\tau (\mu s)$';
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = 25;
ax.XLabel.Position(2) = ax.XLabel.Position(2)+3.5 ; 

ax.YLabel.String = '$\delta$ (MHz)';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = 25;
ax.YLabel.Position(1) = ax.YLabel.Position(1)+0.3; 

cb = colorbar;
caxis([0 1]);
cb.Ticks = 0:0.5:1;

cb.Label.FontSize = 25;


ax.XTick = linspace(3,9,2);
ax.YTick = linspace(-20,20,3);
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