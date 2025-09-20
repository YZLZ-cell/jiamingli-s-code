clear;
clc;

% Parameter setup
t1 = 0.5;
t2 = 5;
B = 100;
w = 4;
p = linspace(pi, 2*pi, 5);
r = linspace(2.5, 4.5, 10);
T = linspace(0.05, 0.2, 6);

% SAP1's F
ave1 = zeros(1, length(T));
final_transfer1 = zeros(length(T), 10);

for i = 1:length(T)
    T_val = T(i);
    d = linspace(0, B/2, 10);
    for j = 1:length(d)
        d_val = d(j); 
        N = round(1000 * (t2 + 2 * t1));
        final_transfer1(i, j) = eve1(N, (0.2-(w-3)*0.02), t1, t2, B * (0.2-(w-3)*0.02) / t2, w, d_val);
    end
    ave1(i) = mean(final_transfer1(i, :));
end

% SAP2's F
ave2 = zeros(1, length(T));
avd2 = zeros(length(T), length(r));
current_avg2 = zeros(length(T), length(r), length(p));
final_transfer2 = zeros(length(T), length(r), length(p), 10);

for i = 1:length(T)
    T_val = T(i);
    d = linspace(0, B/2, 10);
    for m = 1:length(r)
        r_val = r(m);
        for j = 1:length(p)
            p_val = p(j);
            for k = 1:length(d)
                d_val = d(k); 
                N = round(1000 * (t2 + 2 * t1));
                d0_val = B/2 - r_val * 2 * (exp(-t1/T_val) - exp(-(-t1)/T_val)) / (exp(-t1/T_val) + exp(-(-t1)/T_val));
                final_transfer2(i, m, j, k) = eve2(N, T_val, t1, t2,  r_val , B/2 * T_val / t2, w, d_val, d0_val, p_val);
            end
            current_avg2(i, m, j) = mean(final_transfer2(i, m, j, :));
        end
        avd2(i, m) = mean(current_avg2(i, m, :));
    end
    ave2(i) = max(avd2(i, :));
end

% SAP3's F
ave3 = zeros(1, length(T));
avd3 = zeros(length(T), length(r));
current_avg3 = zeros(length(T), length(r), length(p));
final_transfer3 = zeros(length(T), length(r), length(p), 10);

for i = 1:length(T)
    T_val = T(i);
    d = linspace(0, B/2, 10);
    for m = 1:length(r)
        r_val = r(m);
        for j = 1:length(p)
            p_val = p(j);
            for k = 1:length(d)
                d_val = d(k); 
                N = round(1000 * (t2 + 2 * t1));
                d0_val = B/3 - r_val * 2 * (exp(-t1/T_val) - exp(-(-t1)/T_val)) / (exp(-t1/T_val) + exp(-(-t1)/T_val));
                final_transfer3(i, m, j, k) = eve3(N, T_val, t1, t2,  r_val , B/3 * T_val / t2, w, d_val, d0_val, p_val);
            end
            current_avg3(i, m, j) = mean(final_transfer3(i, m, j, :));
        end
        avd3(i, m) = mean(current_avg3(i, m, :));
    end
    ave3(i) = max(avd3(i, :));
end


% Fit curve
coefficients_ave1 = polyfit(T, ave1, 3);
y_fit_ave1 = polyval(coefficients_ave1, T);  

coefficients_ave2 = polyfit(T, ave2, 3);
y_fit_ave2 = polyval(coefficients_ave2, T); 

coefficients_ave3 = polyfit(T, ave3, 3);
y_fit_ave3 = polyval(coefficients_ave3, T);  


figure;
hold on;


h_ave1_scatter = scatter(T, ave1, 75, 'x',...  
    'MarkerEdgeColor', [0.1 0.4 0.7],...  
    'LineWidth', 2.0); 

h_ave1_fit = plot(T, y_fit_ave1, ':',...  
    'Color', [0.2 0.8 0.2],... 
    'LineWidth', 2.5);


h_ave2_scatter = scatter(T, ave2, 75, '+',...
    'MarkerEdgeColor', [0.1 0.5 0.9],...
    'LineWidth', 2.0);

h_ave2_fit = plot(T, y_fit_ave2, '--',...
    'Color', [0.1 0.6 1.0],... 
    'LineWidth', 2.5);


h_ave3_scatter = scatter(T, ave3, 75, 'o',...
    'MarkerFaceColor', [0.9 0.2 0.1],... 
    'MarkerEdgeColor', 'k',...
    'LineWidth', 1.5);

h_ave3_fit = plot(T, y_fit_ave3, '-.',...
    'Color', [1.0 0.3 0.2],...  
    'LineWidth', 2.5);

% Figure setup
ax = gca;
ax.XLabel.String = '$T (\mu s)$ ';
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = 25;
ax.XLabel.Position(2) = ax.XLabel.Position(2)-0.7 ;

ax.YLabel.String = '$F$';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = 25;
ax.YLabel.Position(1) = ax.YLabel.Position(1)+0.02; 
ax.YLabel.Position(2) = ax.YLabel.Position(2)-0.25;

xlim([0.05,0.2]);
ylim([0 1]);

ax.XTick = linspace(0.05, 0.2, 2);  
ax.YTick = linspace(0, 1, 2);           
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 25;
grid off;
box on;

legend([h_ave1_scatter, h_ave2_scatter, h_ave3_scatter], ...
    {'SAP1', 'SAP2', 'SAP3'}, ...
    'Location', 'southeast', ...
    'FontSize', 20, ...
    'Interpreter', 'latex');



hold off;
drawnow;
%% Auxiliary functions
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

function h1 = ham1(t, T, t1, t2, r, w, d)
    a = r * T * log(1 / (exp((t - t1) / T) + exp(-(t - t1) / T))) + ...
        r / T * t2 / 2 * t - ...
        r * T * log(1 / (exp((0 - t1) / T) + exp(-(0 - t1) / T)));
    b = r * T * log(1 / (exp((t - t1 - t2) / T) + exp(-(t - t1 - t2) / T))) - ...
        r / T * t2 / 2 * (t - t2 - t1);
    c = -r / T * (t - t1 - t2 / 2)^2 / 2 + r / T * (-t2 / 2)^2 / 2;
    e = r * T * log(1 / (exp((t1 - t1) / T) + exp(-(t1 - t1) / T))) + ...
        r / T * t2 / 2 * t1 - ...
        r * T * log(1 / (exp((0 - t1) / T) + exp(-(0 - t1) / T)));

    if t < t1
        x = 2 * w / (exp((t - t1) / T) + exp(-(t - t1) / T)) * exp(1i * a );
        y = 2 * w / (exp((t - t1) / T) + exp(-(t - t1) / T)) * exp(-1i * a );
        z = d;
        h1 = x * sig_up() + y * sig_down() + z * sig_z();
    elseif t1 <= t && t < t1 + t2
        x = w * exp(1i * (c + e) );
        y = w * exp(-1i * (c + e));
        z = d;
        h1 = x * sig_up() + y * sig_down() + z * sig_z();
    else
        x = 2 * w / (exp((t - t1 - t2) / T) + exp(-(t - t1 - t2) / T)) * exp(1i * (b + e) );
        y = 2 * w / (exp((t - t1 - t2) / T) + exp(-(t - t1 - t2) / T)) * exp(-1i * (b + e) );
        z = d;
        h1 = x * sig_up() + y * sig_down() + z * sig_z();
    end
end

function result1 = eve1(N, T, t1, t2, r, w, d)
    rho0 = sig_z() + sig_0() / 2;
    rho = rho0;
    result_arr1 = zeros(N, 1);

    for i = 1:N
        dt = 0.001;
        h1 = ham1(i * dt, T, t1, t2, r, w, d);
        rho = rho - 1i * com(h1, rho) * dt - 0.5 * com(h1, com(h1, rho)) * (dt^2);
        result_arr1(i) = real(rho(2, 2));
    end

    result1 = result_arr1(end);
end

function h2 = ham2(t, T, t1, t2, r, r1, w, d, d0, p)
    a = r * T * log(1 / (exp((t - t1) / T) + exp(-(t - t1) / T))) + r1 / T * t2 / 2 * t - ...
        r * T * log(1 / (exp((0 - t1) / T) + exp(-(0 - t1) / T)));
    b = r * T * log(1 / (exp((t - t1 - t2) / T) + exp(-(t - t1 - t2) / T))) - ...
        r1 / T * t2 / 2 * (t - t2 - t1);
    c = -r1 / T * (t - t1 - t2 / 2)^2 / 2 + r1 / T * (-t2 / 2)^2 / 2;
    e = r * T * log(1 / (exp((t1 - t1) / T) + exp(-(t1 - t1) / T))) + ...
        r1 / T * t2 / 2 * t1 - r * T * log(1 / (exp((0 - t1) / T) + exp(-(0 - t1) / T)));

    if t < t1
        x = 2 * w / (exp((t - t1) / T) + exp(-(t - t1) / T)) * (exp(-1i * a + 1i * d0 / 2 * t+1i*p) + 0.95*exp(1i * a - 1i * d0 / 2 * t));
        y = 2 * w / (exp((t - t1) / T) + exp(-(t - t1) / T)) * (exp(1i * a - 1i * d0 / 2 * t-1i*p) + 0.95*exp(-1i * a + 1i * d0 / 2 * t));
        z = d;
    elseif t1 <= t && t < t1 + t2
        x = w * (exp(-1i * (c + e) + 1i * d0 / 2 * t+1i*p) + 0.95*exp(1i * (c + e) - 1i * d0 / 2 * t));
        y = w * (exp(1i * (c + e) - 1i * d0 / 2 * t-1i*p) + 0.95*exp(-1i * (c + e) + 1i * d0 / 2 * t));
        z = d;
    else
        x = 2 * w / (exp((t - t1 - t2) / T) + exp(-(t - t1 - t2) / T)) * (exp(-1i * (b + e) + 1i * d0 / 2 * t+1i*p) + 0.95*exp(1i * (b + e) - 1i * d0 / 2 * t));
        y = 2 * w / (exp((t - t1 - t2) / T) + exp(-(t - t1 - t2) / T)) * (exp(1i * (b + e) - 1i * d0 / 2 * t-1i*p) + 0.95*exp(-1i * (b + e) + 1i * d0 / 2 * t));
        z = d;
    end

    h2 = x * sig_up() + y * sig_down() + z * sig_z();
end


function result2 = eve2(N, T, t1, t2, r, r1, w, d, d0, p)
    rho0 = sig_z() + sig_0() / 2;
    rho = rho0;
    result_arr2 = zeros(N, 1);

    for i = 1:N
        dt = 0.001;
        h2 = ham2(i * dt, T, t1, t2, r, r1, w, d, d0, p);
        rho = rho - 1i * com(h2, rho) * dt - 0.5 * com(h2, com(h2, rho)) * (dt^2);
        result_arr2(i) = real(rho(2, 2));
    end

    result2 = result_arr2(end);
end


function h3 = ham3(t, T, t1, t2, r, r1, w, d, d0, p)
    a = r*T*log(1/(exp((t-t1)/T) + exp(-(t-t1)/T))) + r1/T*t2/2*t - r*T*log(1/(exp((0-t1)/T) + exp(-(0-t1)/T)));
    b = r*T*log(1/(exp((t-t1-t2)/T) + exp(-(t-t1-t2)/T))) - r1/T*t2/2*(t-t2-t1);
    c = -r1/T*(t-t1-t2/2)^2/2 + r1/T*(-t2/2)^2/2;
    e = r*T*log(1/(exp((t1-t1)/T) + exp(-(t1-t1)/T))) + r1/T*t2/2*t1 - r*T*log(1/(exp((0-t1)/T) + exp(-(0-t1)/T)));
    
    if t < t1
        x = 2*w/(exp((t-t1)/T) + exp(-(t-t1)/T)) * (0.95*exp(1i*a+1i*p  ) + 0.9*exp(-1i*a - 1i*d0*t) + exp(-1i*a + 1i*(d0)*t) );
        y = 2*w/(exp((t-t1)/T) + exp(-(t-t1)/T)) * (0.95*exp(-1i*a-1i*p ) + 0.9*exp(1i*a + 1i*d0*t) + exp(1i*a - 1i*(d0)*t) );
        z = d;
    elseif t1 <= t && t < t1 + t2
        x = w * (0.95*exp(1i*(c + e)+1i*p) + 0.9*exp(-1i*(c + e) - 1i*d0*t) + exp(-1i*(c + e) + 1i*d0*t) );
        y = w * (0.95*exp(-1i*(c + e)-1i*p ) + 0.9*exp(1i*(c + e) + 1i*d0*t) + exp(1i*(c + e) - 1i*d0*t) );
        z = d;
    else
        x = 2*w/(exp((t-t1-t2)/T) + exp(-(t-t1-t2)/T)) * (0.95*exp(1i*(b + e)+1i*p) + 0.9*exp(-1i*(b + e) - 1i*d0*t) + exp(-1i*(b + e) + 1i*d0*t) );
        y = 2*w/(exp((t-t1-t2)/T) + exp(-(t-t1-t2)/T)) * (0.95*exp(-1i*(b + e)-1i*p ) + 0.9*exp(1i*(b + e) + 1i*d0*t) + exp(1i*(b + e) - 1i*d0*t) );
        z = d;
    end
    
    h3 = x*sig_up() + y*sig_down() + z*sig_z();
end
function result3 = eve3(N, T, t1, t2, r, r1, w, d, d0, p)
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