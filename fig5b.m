clear;
clc;
%%paramerters
t1 = 0.5;
t2 = linspace(4, 9, 8);
w = 3;
p = linspace(pi, 2*pi, 8);
r = linspace(1.5, 3.5, 6);
T = linspace(0.1, 0.3, 4);


W1_save = zeros(length(t2), 5);
W2_save = zeros(length(t2), 5);
W3_save = zeros(length(t2), 5);

ave1 = zeros(length(t2), 5);
final_transfer1 = zeros(length(t2), 5, 10);

for i = 1:length(t2)
    t2_val = t2(i);
    W = linspace(4*t2_val-15, 4*t2_val+10, 5);
    W1_save(i, :) = W;
    for j = 1:length(W)
        W_val = W(j);
        d = linspace(0, W_val/2, 10);
        for k = 1:length(d)
            d_val = d(k); 
            N = round(1000 * (t2_val + 2 * t1));
            final_transfer1(i, j, k) = eve1(N, 0.2, t1, t2_val, W_val * 0.2 / t2_val, w, d_val);
        end
    end
    ave1(i, j) = mean(final_transfer1(i, j, :));
end


ave2 = zeros(length(t2), 5);
av2 = zeros(length(t2), 5, length(T));
avd2 = zeros(length(t2), 5, length(T), length(r));
current_avg2 = zeros(length(t2), 5, length(T), length(r), length(p));
final_transfer2 = zeros(length(t2), 5, length(T), length(r), length(p), 10);

for i = 1:length(t2)
    t2_val = t2(i);
    W = linspace(8*t2_val-10, 8*t2_val+5, 5);
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
                        N = round(1000 * (t2_val + 2 * t1));
                        d0_val = W_val/2 - r_val * 2 * (exp(-t1/T_val) - exp(-(-t1)/T_val)) / (exp(-t1/T_val) + exp(-(-t1)/T_val));
                        final_transfer2(i, j, n, m, k, l) = eve2(N, T_val, t1, t2_val,  r_val , W_val/2 * T_val / t2_val, w, d_val, d0_val, p_val);
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


ave3 = zeros(length(t2), 5);
av3 = zeros(length(t2), 5, length(T));
avd3 = zeros(length(t2), 5, length(T), length(r));
current_avg3 = zeros(length(t2), 5, length(T), length(r), length(p));
final_transfer3 = zeros(length(t2), 5, length(T), length(r), length(p), 10);

for i = 1:length(t2)
    t2_val = t2(i);
    W = linspace(12*t2_val-10, 12*t2_val+5, 5);
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
                        N = round(1000 * (t2_val + 2 * t1));
                        d0_val = W_val/3 - r_val * 2 * (exp(-t1/T_val) - exp(-(-t1)/T_val)) / (exp(-t1/T_val) + exp(-(-t1)/T_val));
                        final_transfer3(i, j, n, m, k, l) = eve3(N, T_val, t1, t2_val,  r_val , W_val/3 * T_val / t2_val, w, d_val, d0_val, p_val);
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


% figure
figure; 
hold on;
colors = {'k', 'k', 'k'};
markers = {'o', 's', 'd'};
lineStyles = {'-', '--', ':'};

fit_curves = cell(1, 3);
raw_points = cell(1, 3);
legend_handles = gobjects(1, 3); 

for caseIdx = 1:3
    W_points = [];
    t_points = [];
    
    for i = 1:length(t2)
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
        [min_diff, j_idx] = min((0.95 - prob_line));
        if min_diff > 0.01
            j_idx = round(length(prob_line)/2);
        end
        
        W_points(end+1) = W_vals(j_idx);
        t_points(end+1) = t2(i)+2*t1;
    end
    
    [sorted_W, sortIdx] = sort(W_points);
    sorted_t = t_points(sortIdx);
    raw_points{caseIdx} = [sorted_W; sorted_t];
    
    % fit curve
    if length(sorted_W) >= 2
        max_order = min(4, length(sorted_W)-1);
        models = cell(1, max_order);
        gof = zeros(1, max_order);
        for order = 1:max_order
            [p, S] = polyfit(sorted_W, sorted_t, order);
            models{order} = p;
            y_fit = polyval(p, sorted_W);
            y_mean = mean(sorted_t);
            ss_total = sum((sorted_t - y_mean).^2);
            ss_res = sum((sorted_t - y_fit).^2);
            gof(order) = 1 - (ss_res/ss_total);
        end
        [~, best_order] = max(gof);
        best_p = models{best_order};
        fit_x = linspace(min(sorted_W), max(sorted_W), 200);
        fit_y = polyval(best_p, fit_x);
    else
        fit_x = sorted_W;
        fit_y = sorted_t;
    end
    
    fit_curves{caseIdx} = [fit_x; fit_y];
end

%
ax = gca;
ax.XLabel.String = '$W$ (MHz)';
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = 25;
ax.YLabel.String = '$\tau$ ($\mu s$)';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = 25;
xlim([10, 100]);
ylim([5, 10]);  
ax.XTick = linspace(10, 100, 3);
ax.YTick = linspace(5, 10, 2);
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 25;
grid off;

x_fill = linspace(10, 100, 500);
y_curves = zeros(3, length(x_fill));

for i = 1:3
    [sorted_x, sort_idx] = sort(fit_curves{i}(1, :));
    sorted_y = fit_curves{i}(2, sort_idx);
    y_curves(i, :) = interp1(sorted_x, sorted_y, x_fill, 'spline', 'extrap');
    y_curves(i, :) = min(max(y_curves(i, :), 5), 10);
end

boundary_points = cell(1, 4);

boundary_points{1} = [x_fill, fliplr(x_fill); 
                     ones(1, length(x_fill))*5, fliplr(min(y_curves))];

boundary_points{2} = [x_fill, fliplr(x_fill); 
                     min(y_curves), fliplr(median(y_curves))];

boundary_points{3} = [x_fill, fliplr(x_fill); 
                     median(y_curves), fliplr(max(y_curves))];

boundary_points{4} = [x_fill, fliplr(x_fill); 
                     max(y_curves), ones(1, length(x_fill))*10];

fill_colors = {
    [0.8, 0.2, 0.2, 0.7],   
    [0.2, 0.8, 0.2, 0.7],   
    [0.2, 0.2, 0.8, 0.7],   
    [0.8, 0.8, 0.2, 0.7]    
};

for i = 1:4
    poly = polyshape(boundary_points{i}(1, :), boundary_points{i}(2, :));
    poly = simplify(poly);
    plot(poly, 'FaceColor', fill_colors{i}(1:3), ...
               'FaceAlpha', fill_colors{i}(4), ...
               'EdgeColor', 'k', ...
               'LineWidth', 0.5, ...
               'HandleVisibility', 'off');
end

for caseIdx = 1:3
    plot(fit_curves{caseIdx}(1,:), fit_curves{caseIdx}(2,:), ...
         'Color', colors{caseIdx}, ...
         'LineStyle', lineStyles{caseIdx}, ...
         'LineWidth', 3, ...
         'HandleVisibility', 'off'); 
end

for caseIdx = 1:3
    scatter(raw_points{caseIdx}(1,:), raw_points{caseIdx}(2,:), 100, ...
            colors{caseIdx}, markers{caseIdx}, 'filled', ...
            'HandleVisibility', 'off');
end

for caseIdx = 1:3
    legend_handles(caseIdx) = scatter(NaN, NaN, 100, ...
            colors{caseIdx}, markers{caseIdx}, 'filled', ...
            'DisplayName', sprintf('SAP %d', caseIdx));
end

legend(legend_handles, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'southeast');

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