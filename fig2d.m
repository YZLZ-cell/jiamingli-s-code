clear;
clc;
T = 0.5;
t1 = 1;
t2 = 4;
r = 2;
r1 = 2;
w = 4;

%% main
t_values = linspace(0, t2+t1*2, 1000);

y_val1 = r*(exp(-t1/T) - exp(t1/T))/(exp(-t1/T) + exp(t1/T));
y_val2 = -r*(exp(-t1/T) - exp(t1/T))/(exp(-t1/T) + exp(t1/T));

z1_values = zeros(size(t_values));
x1_values = zeros(size(t_values));
z2_values = zeros(size(t_values));
x2_values = zeros(size(t_values));

d0 = r1 / T * t2 - r*2*(exp((-t1)/T) - exp(-(-t1)/T))/(exp((-t1)/T) + exp(-(-t1)/T));

for i = 1:length(t_values)
    t = t_values(i);
    H = ham(t, T, t1, t2, r, r1, w, d0/2);
    z1_values(i) = real(H(1,1));
end

for i = 1:length(t_values)
    t = t_values(i);
    H = ham(t, T, t1, t2, r, r1, w, 0);
    x1_values(i) = real(H(1,2));
end

for i = 1:length(t_values)
    t = t_values(i);
    H = ham1(t, T, t1, t2, r, r1, w, -d0/2);
    z2_values(i) = real(H(1,1));
    x2_values(i) = real(H(1,2));
end

% figure
figure;
hold on;

yyaxis left; 
% Δ₁(t) 
p1 = plot(t_values, z1_values, 'b-', 'LineWidth', 1.5, 'DisplayName', '$\Delta_{1}(t)$');
% Δ₂(t)
p2 = plot(t_values, z2_values, '-', 'LineWidth', 1.5, ...
          'DisplayName', '$\Delta_{2}(t)$', 'Color', [1, 0.5, 0]);

% black dashed lines
hline1 = yline(y_val1, 'k--', 'LineWidth', 1.5, 'DisplayName', '');
hline2 = yline(y_val2, 'k--', 'LineWidth', 1.5, 'DisplayName', '');

yyaxis right; 
% Ω(t)
p3 = plot(t_values, x1_values, 'r-', 'LineWidth', 1.5, 'DisplayName', '$\bar\Omega(t)$');

ax = gca;


yyaxis left;
ax.YLabel.String = '$\delta$ (MHz)';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = 25;
ax.YLabel.Position(1) = ax.YLabel.Position(1)+0.25;
ax.YColor = 'k'; 
ylim([-20, 20]);
yticks(linspace(-20, 20, 3));

yyaxis right;
ax.YLabel.String = '$\bar{\Omega}(t)$ (MHz)';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = 25;
ax.YLabel.Position(1) = ax.YLabel.Position(1)-0.25 ;
ax.YLabel.Position(2) = ax.YLabel.Position(2)-1.6 ;
ax.YColor = 'k'; 
ylim([-20, 20]);
yticks(linspace(-20, 20, 3));

ax.XLabel.String = '$t (\mu s)$';
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = 25;
ax.XLabel.Position(2) = ax.XLabel.Position(2)+3 ; 

xlim([0, 6]);
xticks(linspace(0, 6, 2));

ax.TickLabelInterpreter = 'latex';
ax.FontSize = 25;

ax.Position = [0.13, 0.15, 0.75, 0.75]; 

legend([p1, p2, p3], ...
       'Interpreter', 'latex',...
       'Location', 'southwest',...
       'FontSize', 17.5,...
       'Box', 'on');

grid off;
box on;
%% 
function s = sig_x()
    s = [0 1; 1 0];
end

function s = sig_z()
    s = [1 0; 0 -1];
end

function H = ham(t, T, t1, t2, r, r1, w, d)
    if t < t1
        denom = exp((t-t1)/T) + exp(-(t-t1)/T);
        x = 2*w/denom + d;
        z = r*(exp((t-t1)/T) - exp(-(t-t1)/T))/denom - (r1/T)*(t2/2) + d;
    elseif t < t1 + t2
        x = w + d;
        z = (r1/T)*(t - t1 - t2/2) + d;
    else
        denom = exp((t-t1-t2)/T) + exp(-(t-t1-t2)/T);
        x = 2*w/denom + d;
        z = r*(exp((t-t1-t2)/T) - exp(-(t-t1-t2)/T))/denom + (r1/T)*(t2/2) + d;
    end
    H = x * sig_x() + z * sig_z();
end

function H = ham1(t, T, t1, t2, r, r1, w, d)
    if t < t1
        denom = exp((t-t1)/T) + exp(-(t-t1)/T);
        x = 2*w/denom + d;
        z = -r*(exp((t-t1)/T) - exp(-(t-t1)/T))/denom + (r1/T)*(t2/2) + d;
    elseif t < t1 + t2
        x = w + d;
        z = -(r1/T)*(t - t1 - t2/2) + d;
    else
        denom = exp((t-t1-t2)/T) + exp(-(t-t1-t2)/T);
        x = 2*w/denom + d;
        z = -r*(exp((t-t1-t2)/T) - exp(-(t-t1-t2)/T))/denom - (r1/T)*(t2/2) + d;
    end
    H = x * sig_x() + z * sig_z();
end

