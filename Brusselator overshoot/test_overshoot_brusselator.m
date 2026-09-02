clear; clc; close all;

%% =========================
% 1. USER SETTINGS
% ==========================

% Spatial grid
N = 400;
Lx = 40;
x  = linspace(0, Lx, N)';

% Brusselator / continuation-style parameter
par.sigma = 0.45;

% Overshoot forcing parameters
b0    = 10.0;   % baseline
bmax  = 11.0;   % peak overshoot
Tup   = 40;     % ramp up time
Thold = 80;     % hold time
Tdown = 40;     % ramp down time

% Total simulation time
Tfinal = 300;

% Plot flag for solver
flag = 0;   % keep 0 here; we will make custom plots below

%% =========================
% 2. BUILD FORCING FUNCTION
% ==========================

par.Bfun = @(t) overshoot_B(t, b0, bmax, Tup, Thold, Tdown);

%% =========================
% 3. INITIAL CONDITION
% ==========================

% De Wit equilibrium at baseline b0
a      = 2.5;
u1_ss  = a;
u2_ss  = b0 / a;

% small perturbation around baseline equilibrium
rng(1);
u1 = u1_ss + 0.01 * randn(N,1);
u2 = u2_ss + 0.01 * randn(N,1);

ic0 = [x, u1, u2];

%% =========================
% 4. RUN FORCED SOLVER
% ==========================

[ic_final, S, V] = solve_brusselator_1d_forced(ic0, par, Tfinal, flag);

% Recover full final profiles
u1_final = ic_final(:,2);
u2_final = ic_final(:,3);

%% =========================
% 5. BUILD FORCING CURVE
% ==========================

tt = linspace(0, Tfinal, 1500);
bb = arrayfun(par.Bfun, tt);

%% =========================
% 6. CHOOSE SNAPSHOT TIMES
% ==========================

snapshot_times = [0, 20, 40, 80, 120, 160, 220, 300];
snapshot_inds  = zeros(size(snapshot_times));

for k = 1:length(snapshot_times)
    [~, snapshot_inds(k)] = min(abs(S - snapshot_times(k)));
end

%% =========================
% 7. PLOT FORCING
% ==========================

figure;
plot(tt, bb, 'LineWidth', 1.8);
xlabel('t');
ylabel('b(t)');
title('Overshoot forcing');
box on;

%% =========================
% 8. SPACE-TIME PLOT
% ==========================

figure;
imagesc(1:N, S, V);
axis xy;
xlabel('Spatial index');
ylabel('Time');
title('Space-time plot of u_2(x,t)');
colorbar;
box on;

%% =========================
% 9. SNAPSHOT PROFILES
% ==========================

figure;
hold on;
for k = 1:length(snapshot_inds)
    plot(x, V(snapshot_inds(k),:), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('t = %.1f', S(snapshot_inds(k))));
end
xlabel('x');
ylabel('u_2(x,t)');
title('Snapshots of u_2 during overshoot experiment');
legend('Location','best');
box on;
hold off;

%% =========================
% 10. FINAL PROFILE
% ==========================

figure;
plot(x, u2_final, 'LineWidth', 1.8);
xlabel('x');
ylabel('u_2(x,T)');
title(sprintf('Final profile at T = %.1f', Tfinal));
box on;

%% =========================
% 11. SAVE OUTPUT
% ==========================

run_data.ic0            = ic0;
run_data.ic_final       = ic_final;
run_data.x              = x;
run_data.S              = S;
run_data.V              = V;
run_data.b0             = b0;
run_data.bmax           = bmax;
run_data.Tup            = Tup;
run_data.Thold          = Thold;
run_data.Tdown          = Tdown;
run_data.Tfinal         = Tfinal;
run_data.sigma          = par.sigma;
run_data.snapshot_times = snapshot_times;
run_data.snapshot_inds  = snapshot_inds;

save('overshoot_test_run.mat', 'run_data');

disp('Done. Saved results to overshoot_test_run.mat');