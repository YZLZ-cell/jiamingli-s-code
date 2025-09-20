clear;
clc;

%%paramerters

t1 = 0.5;
t2 = 4;
w = linspace(1, 5, 8);
p = linspace(pi, 2*pi, 10);
T = linspace(0.1, 0.3, 3);

W1_save = zeros(length(w), 5);
W2_save = zeros(length(w), 5);
W3_save = zeros(length(w), 5);

ave1 = zeros(length(w), 5);
final_transfer1 = zeros(length(w), 5, 10);

for i = 1:length(w)
    w_val = w(i);
    W = linspace(3*4*w_val*w_val/9, 3*4*w_val*w_val/9+5, 5);
    W1_save(i, :) = W;
    for j = 1:length(W)
        W_val = W(j);
        d = linspace(0, W_val/2, 10);
        for k = 1:length(d)
            d_val = d(k); 
            N = round(1000 * (t2 + 2 * t1));
            final_transfer1(i, j, k) = eve1(N, 0.3, t1, t2, W_val * 0.3 / t2, w_val, d_val);
        end
    end
    ave1(i, j) = mean(final_transfer1(i, j, :));
end

ave2 = zeros(length(w), 5);
av2 = zeros(length(w), 5, length(T));
avd2 = zeros(length(w), 5, length(T), 5);
current_avg2 = zeros(length(w), 5, length(T), 5, length(p));
final_transfer2 = zeros(length(w), 5, length(T), 5, length(p), length(d));

for i = 1:length(w)
    w_val = w(i);
    r = linspace(w_val-2, w_val, 5);
    W = linspace(6*4*w_val*w_val/9+10, 6*4*w_val*w_val/9+15, 5);
    W2_save(i, :) = W;
    for j = 1:length(W)
        W_val = W(j);
        d = linspace(0, W_val/2, 10);
        for n = 1:length(T)
            T_val = T(n);
            for m = 1:length(r)
                r_val = r(m);
                for k = 1:length(p)
                    p_val = p(k);
                    for l = 1:length(d)
                        d_val = d(l); 
                        N = round(1000 * (t2 + 2 * t1));
                        d0_val = W_val/2.1 - r_val * 2 * (exp(-t1/T_val) - exp(-(-t1)/T_val)) / (exp(-t1/T_val) + exp(-(-t1)/T_val));
                        final_transfer2(i, j, n, m, k, l) = eve2(N, T_val, t1, t2,  r_val , W_val/2.1 * T_val / t2, w_val, d_val, d0_val, p_val);
                    end
                    current_avg2(i, j, n, m,  k) = mean(final_transfer2(i, j, n, m, k, :));
                end
                avd2(i, j, n, m) = mean(current_avg2(i, j, n, m, :));
            end
            av2(i, j, n) = max(avd2(i, j, n, :));
        end
        ave2(i, j) = max(av2(i, j, :));
    end
end

ave3 = zeros(length(w), 5);
av3 = zeros(length(w), 5, length(T));
avd3 = zeros(length(w), 5, length(T), 5);
current_avg3 = zeros(length(w), 5, length(T), 5, length(p));
final_transfer3 = zeros(length(w), 5, length(T), 5, length(p), length(d));

for i = 1:length(w)
    w_val = w(i);
    r = linspace(w_val-2, w_val, 5);
    W = linspace(9*4*w_val*w_val/9+15, 9*4*w_val*w_val/9+20, 5);
    W3_save(i, :) = W;
    for j = 1:length(W)
        W_val = W(j);
        d = linspace(0, W_val/2, 10);
        for n = 1:length(T)
            T_val = T(n);
            for m = 1:length(r)
                r_val = r(m);
                for k = 1:length(p)
                    p_val = p(k);
                    for l = 1:length(d)
                        d_val = d(l); 
                        N = round(1000 * (t2 + 2 * t1));
                        d0_val = W_val/3.1 - r_val * 2 * (exp(-t1/T_val) - exp(-(-t1)/T_val)) / (exp(-t1/T_val) + exp(-(-t1)/T_val));
                        final_transfer3(i, j, n, m, k, l) = eve3(N, T_val, t1, t2,  r_val , W_val/3.1 * T_val / t2, w_val, d_val, d0_val, p_val);
                    end
                    current_avg3(i, j, n, m, k) = mean(final_transfer3(i, j, n, m, k, :));
                end
                avd3(i, j, n, m) = mean(current_avg3(i, j, n, m, :));
            end
            av3(i, j, n) = max(avd3(i, j, n, :));
        end
        ave3(i, j) = max(av3(i, j, :));
    end
end   

figure; hold on;
colors = {'r', 'g', 'b'};
markers = {'o', 's', 'd'};
lineStyles = {'-', '--', ':'};

fit_curves = cell(1, 3);

for caseIdx = 1:3
    W_points = [];
    w_points = [];
    
    for i = 1:length(w)
        switch caseIdx
            case 1
                prob_line = ave1(i, :);
                W_vals = W1_save(i, :);
            case 2
                prob_line = ave2(i, :);
                W_vals = W2_save(i, :);
            case 3
                prob_line = ave3(i, :);
                W_vals = W3_save(i, :);
        end
        
        % find the points nearest to 0.95
        [~, j_idx] = min(abs(prob_line - 0.95));
      
        W_points(end+1) = W_vals(j_idx);
        w_points(end+1) = w(i);
    end
    
    [sorted_W, sortIdx] = sort(W_points);
    sorted_w = w_points(sortIdx);
    
    % fit curve
    if length(sorted_W) > 2
        p = polyfit(sorted_W, sorted_w, 2);
        fit_x = linspace(min(sorted_W), max(sorted_W), 100);
        fit_y = polyval(p, fit_x);
    else
        fit_x = linspace(min(sorted_W), max(sorted_W), 100);
        fit_y = interp1(sorted_W, sorted_w, fit_x, 'linear', 'extrap');
    end
    fit_curves{caseIdx} = [fit_x; fit_y];
    plot(fit_x, fit_y, ...
        'Color', colors{caseIdx}, ...
        'LineStyle', lineStyles{caseIdx}, ...
        'LineWidth', 3, ...
        'DisplayName', sprintf('RPA %d', caseIdx));
    scatter(sorted_W, sorted_w, 100, colors{caseIdx}, markers{caseIdx}, 'filled', ...
        'HandleVisibility', 'off');
end

% fill colors

x1 = fit_curves{1}(1,:);
y1 = fit_curves{1}(2,:);
x2 = fit_curves{2}(1,:);
y2 = fit_curves{2}(2,:);
x3 = fit_curves{3}(1,:);
y3 = fit_curves{3}(2,:);

y_ref = linspace(1, 5, 500);
x1_interp = interp1(y1, x1, y_ref, 'linear', 'extrap');
x2_interp = interp1(y2, x2, y_ref, 'linear', 'extrap');
x3_interp = interp1(y3, x3, y_ref, 'linear', 'extrap');

x_min = min([x1_interp, x2_interp, x3_interp]);
x_max = max([x1_interp, x2_interp, x3_interp]);
x1_interp = max(x1_interp, x_min);
x1_interp = min(x1_interp, x_max);
x2_interp = max(x2_interp, x_min);
x2_interp = min(x2_interp, x_max);
x3_interp = max(x3_interp, x_min);
x3_interp = min(x3_interp, x_max);

x1_interp(isnan(x1_interp)) = x_min;
x2_interp(isnan(x2_interp)) = x_min;
x3_interp(isnan(x3_interp)) = x_min;

colors_fill = {
    [0.9, 0.9, 1.0, 0.3]   
    [0.9, 1.0, 0.9, 0.3]  
    [1.0, 0.9, 0.9, 0.3]   
    [0.95, 0.95, 0.8, 0.3] 
};

X_region1 = [x_min*ones(size(y_ref)), fliplr(x1_interp)];
Y_region1 = [y_ref, fliplr(y_ref)];
fill(X_region1, Y_region1, colors_fill{1}, 'EdgeColor', 'none', 'HandleVisibility', 'off');

X_region2 = [x1_interp, fliplr(x2_interp)];
Y_region2 = [y_ref, fliplr(y_ref)];
fill(X_region2, Y_region2, colors_fill{2}, 'EdgeColor', 'none', 'HandleVisibility', 'off');

X_region3 = [x2_interp, fliplr(x3_interp)];
Y_region3 = [y_ref, fliplr(y_ref)];
fill(X_region3, Y_region3, colors_fill{3}, 'EdgeColor', 'none', 'HandleVisibility', 'off');

X_region4 = [x3_interp, x_max*ones(size(y_ref))];
Y_region4 = [y_ref, fliplr(y_ref)];
fill(X_region4, Y_region4, colors_fill{4}, 'EdgeColor', 'none', 'HandleVisibility', 'off');

for caseIdx = 1:3
    fit_x = fit_curves{caseIdx}(1,:);
    fit_y = fit_curves{caseIdx}(2,:);
    plot(fit_x, fit_y, ...
        'Color', colors{caseIdx}, ...
        'LineStyle', lineStyles{caseIdx}, ...
        'LineWidth', 3);
    
    switch caseIdx
        case 1
            for i = 1:length(w)
                [~, j_idx] = min(abs(ave1(i, :) - 0.95));
                scatter(W1_save(i, j_idx), w(i), 100, colors{caseIdx}, markers{caseIdx}, 'filled');
            end
        case 2
            for i = 1:length(w)
                [~, j_idx] = min(abs(ave2(i, :) - 0.95));
                scatter(W2_save(i, j_idx), w(i), 100, colors{caseIdx}, markers{caseIdx}, 'filled');
            end
        case 3
            for i = 1:length(w)
                [~, j_idx] = min(abs(ave3(i, :) - 0.95));
                scatter(W3_save(i, j_idx), w(i), 100, colors{caseIdx}, markers{caseIdx}, 'filled');
            end
    end
end

% 
ax = gca;
ax.XLabel.String = '$W$ (MHz)';
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = 25;
ax.YLabel.String = '$\Omega$ (MHz)';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = 25;

xlim([10, 100]);
ylim([1, 5]);

ax.XTick = linspace(10, 100, 3);
ax.YTick = linspace(1, 5, 2);
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 25;
grid off;

legend('show', 'Interpreter', 'latex', 'FontSize', 15, 'Location', 'best');

box on;
set(gcf, 'Color', 'w');
set(gca, 'LineWidth', 2);


hold off;
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
        x = 2 * w / (exp((t - t1) / T) + exp(-(t - t1) / T)) * (exp(-1i * a + 1i * d0 / 2 * t+1i*p) + exp(1i * a - 1i * d0 / 2 * t));
        y = 2 * w / (exp((t - t1) / T) + exp(-(t - t1) / T)) * (exp(1i * a - 1i * d0 / 2 * t-1i*p) + exp(-1i * a + 1i * d0 / 2 * t));
        z = d;
    elseif t1 <= t && t < t1 + t2
        x = w * (exp(-1i * (c + e) + 1i * d0 / 2 * t+1i*p) + exp(1i * (c + e) - 1i * d0 / 2 * t));
        y = w * (exp(1i * (c + e) - 1i * d0 / 2 * t-1i*p) + exp(-1i * (c + e) + 1i * d0 / 2 * t));
        z = d;
    else
        x = 2 * w / (exp((t - t1 - t2) / T) + exp(-(t - t1 - t2) / T)) * (exp(-1i * (b + e) + 1i * d0 / 2 * t+1i*p) + exp(1i * (b + e) - 1i * d0 / 2 * t));
        y = 2 * w / (exp((t - t1 - t2) / T) + exp(-(t - t1 - t2) / T)) * (exp(1i * (b + e) - 1i * d0 / 2 * t-1i*p) + exp(-1i * (b + e) + 1i * d0 / 2 * t));
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
        x = 2*w/(exp((t-t1)/T) + exp(-(t-t1)/T)) * (exp(1i*a  ) + exp(-1i*a - 1i*d0*t+1i*p) + exp(-1i*a + 1i*(d0)*t) );
        y = 2*w/(exp((t-t1)/T) + exp(-(t-t1)/T)) * (exp(-1i*a ) + exp(1i*a + 1i*d0*t-1i*p) + exp(1i*a - 1i*(d0)*t) );
        z = d;
    elseif t1 <= t && t < t1 + t2
        x = w * (exp(1i*(c + e)) + exp(-1i*(c + e) - 1i*d0*t+1i*p) + exp(-1i*(c + e) + 1i*d0*t) );
        y = w * (exp(-1i*(c + e) ) + exp(1i*(c + e) + 1i*d0*t-1i*p) + exp(1i*(c + e) - 1i*d0*t) );
        z = d;
    else
        x = 2*w/(exp((t-t1-t2)/T) + exp(-(t-t1-t2)/T)) * (exp(1i*(b + e)) + exp(-1i*(b + e) - 1i*d0*t+1i*p) + exp(-1i*(b + e) + 1i*d0*t) );
        y = 2*w/(exp((t-t1-t2)/T) + exp(-(t-t1-t2)/T)) * (exp(-1i*(b + e) ) + exp(1i*(b + e) + 1i*d0*t-1i*p) + exp(1i*(b + e) - 1i*d0*t) );
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