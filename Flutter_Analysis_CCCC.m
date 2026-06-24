%% Flutter Analysis for PZT/0/90/90/0/PZT Laminate with CCCC Boundary
% Based on Tawfik's panel flutter formulation (Section 9.1-9.2)
% Plotting style from Roberts (2022) thesis

%% ========================================================================
% SECTION 1: PLATE GEOMETRY AND MATERIAL PROPERTIES
% ========================================================================
a = 0.200;           % Length (m) - flow direction
b = 0.300;           % Width (m) - span direction
area = a * b;

fprintf('=== PLATE DIMENSIONS ===\n');
fprintf('Length a (flow direction) = %.4f m\n', a);
fprintf('Width b = %.4f m\n', b);
fprintf('Aspect ratio a/b = %.2f\n', a/b);
fprintf('\n');

%% Material definitions
pzt = struct('type', 'piezo', ...
             'E', 66e9, ...
             'nu', 0.31, ...
             'rho', 7800, ...
             'thickness', 0.00015, ...
             'd31',-171e-12,...
             'e31', -5.4, ...
             'eps33', 1.5e-8);

comp0 = struct('type', 'comp', ...
               'E1', 140e9, ...
               'E2', 10e9, ...
               'G12', 5e9, ...
               'nu12', 0.3, ...
               'rho', 1600, ...
               'thickness', 0.00025, ...
               'angle', 0);

comp90 = comp0;
comp90.angle = 90;

layers = {pzt, comp0, comp90, comp90, comp0, pzt};

%% Calculate laminate properties
[A_mat, B, D, As, I0, I2] = calculate_laminate_properties(layers);

D11 = D(1,1);
D22 = D(2,2);
D12 = D(1,2);
D66 = D(3,3);
I0_mass = I0;
h_total = sum(cellfun(@(x) x.thickness, layers));

fprintf('=== Laminate Properties ===\n');
fprintf('D11 (bending stiffness) = %.2f N·m\n', D11);
fprintf('D22 = %.2f N·m\n', D22);
fprintf('I0 = %.4f kg/m²\n', I0_mass);
fprintf('Total thickness h = %.4f m (%.1f mm)\n', h_total, h_total*1000);
fprintf('\n');

%% ========================================================================
% SECTION 2: CORRECT CCCC BEAM EIGENVALUES (All Edges Clamped)
% ========================================================================
% x-direction: Clamped-Clamped (C-C)
kL_CC = [4.73004074486270; 7.85320462409584; 10.9956078380017];
k_x1 = kL_CC(1)/a;
k_x2 = kL_CC(2)/a;

sigma_x1 = (cosh(kL_CC(1)) - cos(kL_CC(1))) / (sinh(kL_CC(1)) - sin(kL_CC(1)));
sigma_x2 = (cosh(kL_CC(2)) - cos(kL_CC(2))) / (sinh(kL_CC(2)) - sin(kL_CC(2)));

phi1 = @(x) cosh(k_x1*x) - cos(k_x1*x) - sigma_x1*(sinh(k_x1*x) - sin(k_x1*x));
phi2 = @(x) cosh(k_x2*x) - cos(k_x2*x) - sigma_x2*(sinh(k_x2*x) - sin(k_x2*x));

% y-direction: Clamped-Clamped (C-C)
kL_CC_y = [4.73004074486270; 7.85320462409584];
k_y1 = kL_CC_y(1)/b;

sigma_y1 = (cosh(kL_CC_y(1)) - cos(kL_CC_y(1))) / (sinh(kL_CC_y(1)) - sin(kL_CC_y(1)));
psi1 = @(y) cosh(k_y1*y) - cos(k_y1*y) - sigma_y1*(sinh(k_y1*y) - sin(k_y1*y));

% Mass normalization
norm1 = sqrt(integral2(@(x,y) I0_mass * (phi1(x).*psi1(y)).^2, 0, a, 0, b));
norm2 = sqrt(integral2(@(x,y) I0_mass * (phi2(x).*psi1(y)).^2, 0, a, 0, b));

phi1_n = @(x) phi1(x)/sqrt(norm1/I0_mass);
phi2_n = @(x) phi2(x)/sqrt(norm2/I0_mass);
psi1_n = @(y) psi1(y)/sqrt(b);

mode1 = @(x,y) phi1_n(x) .* psi1_n(y);
mode2 = @(x,y) phi2_n(x) .* psi1_n(y);

% VERIFY boundary conditions
fprintf('\n=== CCCC BOUNDARY CONDITION VERIFICATION ===\n');
fprintf('phi1(0) = %.6f (should be 0)\n', phi1(0));
fprintf('phi1(a) = %.6f (should be 0 for CCCC!)\n', phi1(a));
fprintf('psi1(0) = %.6f (should be 0)\n', psi1(0));
fprintf('psi1(b) = %.6f (should be 0)\n', psi1(b));

%% ========================================================================
% SECTION 3: NUMERICAL INTEGRATION
% ========================================================================
n_quad = 120;
[x_quad, w_x] = gauss_legendre(n_quad, 0, a);
[y_quad, w_y] = gauss_legendre(n_quad, 0, b);
[W_X, W_Y] = meshgrid(w_x, w_y);
W_2D = W_X .* W_Y;

%% ========================================================================
% SECTION 4: AERODYNAMIC INTEGRALS
% ========================================================================
fprintf('Computing aerodynamic integrals...\n');

dx_small = 1e-7;
dphi1 = @(x) (phi1_n(x+dx_small) - phi1_n(x-dx_small))/(2*dx_small);
dphi2 = @(x) (phi2_n(x+dx_small) - phi2_n(x-dx_small))/(2*dx_small);

A_aero = zeros(2,2);
B_aero = zeros(2,2);

for i = 1:n_quad
    for j = 1:n_quad
        x = x_quad(i);
        y = y_quad(j);
        w = W_2D(j,i);
        
        A_aero(1,1) = A_aero(1,1) + mode1(x,y) * dphi1(x) * psi1_n(y) * w;
        A_aero(1,2) = A_aero(1,2) + mode1(x,y) * dphi2(x) * psi1_n(y) * w;
        A_aero(2,1) = A_aero(2,1) + mode2(x,y) * dphi1(x) * psi1_n(y) * w;
        A_aero(2,2) = A_aero(2,2) + mode2(x,y) * dphi2(x) * psi1_n(y) * w;
        
        B_aero(1,1) = B_aero(1,1) + mode1(x,y)^2 * w;
        B_aero(1,2) = B_aero(1,2) + mode1(x,y) * mode2(x,y) * w;
        B_aero(2,1) = B_aero(2,1) + mode2(x,y) * mode1(x,y) * w;
        B_aero(2,2) = B_aero(2,2) + mode2(x,y)^2 * w;
    end
end

%% ========================================================================
% SECTION 5: STRUCTURAL MATRICES
% ========================================================================
fprintf('Computing structural matrices...\n');

M_modal = I0_mass * B_aero;

d2phi1 = @(x) (phi1_n(x+dx_small) - 2*phi1_n(x) + phi1_n(x-dx_small))/(dx_small^2);
d2phi2 = @(x) (phi2_n(x+dx_small) - 2*phi2_n(x) + phi2_n(x-dx_small))/(dx_small^2);
d2psi1 = @(y) (psi1_n(y+dx_small) - 2*psi1_n(y) + psi1_n(y-dx_small))/(dx_small^2);
dpsi1 = @(y) (psi1_n(y+dx_small) - psi1_n(y-dx_small))/(2*dx_small);

K_modal = zeros(2,2);

for i = 1:n_quad
    for j = 1:n_quad
        x = x_quad(i);
        y = y_quad(j);
        w = W_2D(j,i);
        
        w1_xx = d2phi1(x) * psi1_n(y);
        w1_yy = phi1_n(x) * d2psi1(y);
        w1_xy = dphi1(x) * dpsi1(y);
        
        w2_xx = d2phi2(x) * psi1_n(y);
        w2_yy = phi2_n(x) * d2psi1(y);
        w2_xy = dphi2(x) * dpsi1(y);
        
        K_modal(1,1) = K_modal(1,1) + (D11*w1_xx^2 + D22*w1_yy^2 + 2*D12*w1_xx*w1_yy + 4*D66*w1_xy^2) * w;
        K_modal(2,2) = K_modal(2,2) + (D11*w2_xx^2 + D22*w2_yy^2 + 2*D12*w2_xx*w2_yy + 4*D66*w2_xy^2) * w;
        K_modal(1,2) = K_modal(1,2) + (D11*w1_xx*w2_xx + D22*w1_yy*w2_yy + D12*(w1_xx*w2_yy + w1_yy*w2_xx) + 4*D66*w1_xy*w2_xy) * w;
    end
end
K_modal(2,1) = K_modal(1,2);

% Solve eigenvalue problem for natural frequencies
[V_mode, D_mode] = eig(K_modal, M_modal);
omega_n = sqrt(diag(D_mode));
[omega_n, idx] = sort(omega_n);
freq_n = omega_n/(2*pi);

fprintf('Natural frequencies (no flow):\n');
fprintf('  f1 (bending) = %.1f Hz\n', freq_n(1));
fprintf('  f2 (torsion) = %.1f Hz\n', freq_n(2));
fprintf('  f2/f1 = %.2f\n', freq_n(2)/freq_n(1));
fprintf('\n');

%% ========================================================================
% SECTION 6: FLUTTER ANALYSIS (Tawfik Method)
% ========================================================================
% Flow conditions
M_inf = 1.5;
beta_flow = sqrt(M_inf^2 - 1);
rho_air = 1.225;
speed_of_sound = 340;

% Reference frequency
omega_0 = omega_n(1);
L_panel = a;
rho_panel = I0_mass / h_total;

% Thermal load (ΔT = 0)
delta_T = 0;
alpha_thermal = 22.5e-6;
E_avg = (comp0.E1 + comp0.E2)/2;
nu_avg = comp0.nu12;
N_T = alpha_thermal * E_avg * h_total * delta_T / (1 - nu_avg);
K_T = N_T * B_aero;  % Thermal stiffness matrix

% λ range
lambda_min = 0;
lambda_max = 1000;
n_lambda = 400;
lambda_values = linspace(lambda_min, lambda_max, n_lambda);

fprintf('=== FLUTTER ANALYSIS (Tawfik Method) ===\n');
fprintf('Mach = %.2f, β = %.3f\n', M_inf, beta_flow);
fprintf('Panel length L = %.4f m\n', L_panel);
fprintf('λ range: %.0f to %.0f (%d points)\n', lambda_min, lambda_max, n_lambda);
fprintf('\n');

% Storage for YOUR preferred plotting style
freq_Hz = zeros(2, n_lambda);
damping_ratio = zeros(2, n_lambda);
U_velocity = zeros(1, n_lambda);

flutter_detected = false;
lambda_cr = 0;
U_flutter = 0;
freq_flutter = 0;
Mach_flutter = 0;

for i = 1:n_lambda
    lambda = lambda_values(i);
    
    % Velocity from λ
    q_dyn = (lambda * beta_flow * D11) / (2 * L_panel^3);
    if q_dyn > 0
        U = sqrt(2 * q_dyn / rho_air);
    else
        U = 0;
    end
    U_velocity(i) = U;
    
    % Aerodynamic damping
    if U > 0 && M_inf^2 > 2
        g_a = (rho_air * U * (M_inf^2 - 2)) / (rho_panel * h_total * omega_0 * beta_flow^3);
    else
        g_a = 0;
    end
    
    % Matrices
    K_aero = lambda * (D11 / L_panel^3) * A_aero;
    G_aero = g_a * omega_0 * M_modal;
    K_total = K_modal - K_T + K_aero;
    
    % State-space eigenvalue problem
    A_state = [zeros(2), eye(2); -M_modal\K_total, -M_modal\G_aero];
    [~, D_state] = eig(A_state);
    eigenvalues = diag(D_state);
    
    % Store real parts
    real_eigenvalues(1, i) = max(real(eigenvalues));
    real_eigenvalues(2, i) = min(real(eigenvalues));

    % Extract oscillatory modes
    eig_osc = eigenvalues(imag(eigenvalues) > 1e-6);
    [~, sort_idx] = sort(imag(eig_osc));
    eig_osc = eig_osc(sort_idx);
    
    for j = 1:min(2, length(eig_osc))
        freq_Hz(j, i) = imag(eig_osc(j)) / (2*pi);
        % Damping ratio = -real/|ω|
        damping_ratio(j, i) = -real(eig_osc(j)) / abs(eig_osc(j));
    end
    
    % Flutter detection (when damping crosses zero)
    if i > 1 && ~flutter_detected && damping_ratio(2, i) > 0 && damping_ratio(2, i-1) < 0
        flutter_detected = true;
        
        % Linear interpolation for exact crossing
        l1 = lambda_values(i-1);
        l2 = lambda_values(i);
        d1 = damping_ratio(2, i-1);
        d2 = damping_ratio(2, i);
        lambda_cr = l1 - d1 * (l2 - l1) / (d2 - d1);
        
        % Interpolate frequency and velocity
        f1 = freq_Hz(2, i-1);
        f2 = freq_Hz(2, i);
        freq_flutter = f1 + (f2 - f1) * (lambda_cr - l1) / (l2 - l1);
        
        U_flutter = sqrt((lambda_cr * beta_flow * D11) / (rho_air * a^3));
        Mach_flutter = U_flutter / speed_of_sound;
    end
    
    if mod(i, 50) == 0
        fprintf('  Progress: %.0f%% (U = %.1f m/s)\n', i/n_lambda*100, U);
    end
end

% Display results
if flutter_detected
    fprintf('\n========================================\n');
    fprintf('      FLUTTER ANALYSIS RESULTS\n');
    fprintf('========================================\n');
    fprintf('Critical λ_cr: %.1f\n', lambda_cr);
    fprintf('Flutter velocity: %.1f m/s (Mach %.2f)\n', U_flutter, Mach_flutter);
    fprintf('Flutter frequency: %.1f Hz\n', freq_flutter);
    fprintf('========================================\n');
else
    fprintf('\n========================================\n');
    fprintf('      FLUTTER ANALYSIS RESULTS\n');
    fprintf('========================================\n');
    fprintf('No flutter detected within λ range [%.0f, %.0f]\n', lambda_min, lambda_max);
    fprintf('========================================\n');
end

%% ========================================================================
% SECTION 7: PLOTTING (YOUR PREFERRED STYLE)
% ========================================================================
figure('Position', [100, 100, 1200, 500]);

% Subplot 1: Frequency Coalescence
subplot(1,2,1);
plot(lambda_values, freq_Hz(1,:)/1000, 'b-', 'LineWidth', 2); hold on;
plot(lambda_values, freq_Hz(2,:)/1000, 'r-', 'LineWidth', 2);
if flutter_detected
    plot(lambda_cr, freq_flutter/1000, 'ko', 'MarkerSize', 10, 'LineWidth', 2);
    xline(lambda_cr, 'k--');
end
xlabel('\lambda', 'FontSize', 12);
ylabel('Frequency (kHz)', 'FontSize', 12);
title('Frequency Coalescence', 'FontSize', 14);
legend('Mode 1', 'Mode 2', 'Location', 'best');
grid on;

% Subplot 2: Damping Evolution
subplot(1,2,2);
plot(lambda_values, damping_ratio(1,:), 'b-', 'LineWidth', 2); hold on;
plot(lambda_values, damping_ratio(2,:), 'r-', 'LineWidth', 2);
if flutter_detected
    plot(lambda_cr, 0, 'ko', 'MarkerSize', 10, 'LineWidth', 2);
    xline(lambda_cr, 'k--');
end
yline(0, 'k--');
xlabel('\lambda', 'FontSize', 12);
ylabel('Damping Ratio', 'FontSize', 12);
title('Damping Evolution', 'FontSize', 14);
legend('Mode 1', 'Mode 2', 'Location', 'best');
grid on;

if flutter_detected
    sgtitle(sprintf('Flutter: λ_cr = %.1f, f = %.0f Hz (Mach %.2f)', lambda_cr, freq_flutter, Mach_flutter), 'FontSize', 14);
else
    sgtitle('No Flutter Detected', 'FontSize', 14);
end

% Save figure
saveas(gcf, 'flutter_results.png');
saveas(gcf, 'flutter_results.fig');
fprintf('\n✓ Figure saved as: flutter_results.png\n');

%% ========================================================================
% SECTION 7: IMPROVED MODE SHAPES PLOTTING FOR CCCC (All Edges Clamped)
% ========================================================================
figure('Position', [100, 100, 1600, 700]);

% Create fine mesh for plotting
n_plot = 80;
x_plot = linspace(0, a, n_plot);
y_plot = linspace(0, b, n_plot);
[X_plot, Y_plot] = meshgrid(x_plot, y_plot);

% Calculate individual mode shapes
Z_mode1 = zeros(n_plot);
Z_mode2 = zeros(n_plot);
for i = 1:n_plot
    for j = 1:n_plot
        Z_mode1(i,j) = mode1(x_plot(j), y_plot(i));
        Z_mode2(i,j) = mode2(x_plot(j), y_plot(i));
    end
end

% Normalize for visualization
Z_mode1 = Z_mode1 / max(abs(Z_mode1(:)));
Z_mode2 = Z_mode2 / max(abs(Z_mode2(:)));

% ========================================================================
% ROW 1: INDIVIDUAL MODES BEFORE FLUTTER
% ========================================================================

% Mode 1: First Bending Mode
subplot(2,3,1);
surf(X_plot, Y_plot, Z_mode1, 'EdgeColor', 'none');
colormap(jet);
colorbar;
xlabel('Length x (m)', 'FontSize', 10);
ylabel('Width y (m)', 'FontSize', 10);
zlabel('Amplitude', 'FontSize', 10);
title('Mode 1: First Bending Mode (Before Flutter)', 'FontSize', 11, 'FontWeight', 'bold');
view(45, 30);
grid on;
lighting gouraud;
light('Position', [1 1 1]);
hold on;
% Mark all clamped edges (CCCC - all edges clamped)
plot3([0 a a 0 0], [0 0 b b 0], [0 0 0 0 0], 'r-', 'LineWidth', 2);
text(a/2, -0.03, 1.2, 'CLAMPED EDGE', 'FontSize', 8, 'Color', 'r', 'HorizontalAlignment', 'center');
text(a/2, b+0.02, 1.2, 'CLAMPED EDGE', 'FontSize', 8, 'Color', 'r', 'HorizontalAlignment', 'center');
text(-0.05, b/2, 1.2, 'CLAMPED', 'Rotation', 90, 'FontSize', 8, 'Color', 'r');
text(a+0.03, b/2, 1.2, 'CLAMPED', 'Rotation', -90, 'FontSize', 8, 'Color', 'r');

% Mode 2: Torsion/Shear Mode
subplot(2,3,2);
surf(X_plot, Y_plot, Z_mode2, 'EdgeColor', 'none');
colormap(jet);
colorbar;
xlabel('Length x (m)', 'FontSize', 10);
ylabel('Width y (m)', 'FontSize', 10);
zlabel('Amplitude', 'FontSize', 10);
title('Mode 2: Torsion/Shear Mode (Before Flutter)', 'FontSize', 11, 'FontWeight', 'bold');
view(45, 30);
grid on;
lighting gouraud;
hold on;
plot3([0 a a 0 0], [0 0 b b 0], [0 0 0 0 0], 'r-', 'LineWidth', 2);

% ========================================================================
% ROW 1, COL 3: COALESCED MODE AT FLUTTER (50/50 MIX)
% ========================================================================
subplot(2,3,3);
% At flutter, modes coalesce - equal contribution from both modes
Z_coalesced = 0.5 * Z_mode1 + 0.5 * Z_mode2;
Z_coalesced = Z_coalesced / max(abs(Z_coalesced(:)));

surf(X_plot, Y_plot, Z_coalesced, 'EdgeColor', 'none');
colormap(jet);
colorbar;
xlabel('Length x (m)', 'FontSize', 10);
ylabel('Width y (m)', 'FontSize', 10);
zlabel('Amplitude', 'FontSize', 10);
if flutter_detected
    title(sprintf('COALESCED MODE AT FLUTTER (λ_{cr}=%.0f, M=%.2f)', lambda_cr, Mach_flutter), ...
          'FontSize', 11, 'FontWeight', 'bold', 'Color', 'r');
else
    title('COALESCED MODE AT FLUTTER (Theoretical)', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'r');
end
view(45, 30);
grid on;
lighting gouraud;
hold on;
plot3([0 a a 0 0], [0 0 b b 0], [0 0 0 0 0], 'k-', 'LineWidth', 1.5);

% ========================================================================
% ROW 2, COL 1: MODE CONTRIBUTION EVOLUTION (COALESCENCE)
% ========================================================================
subplot(2,3,4);
% Show how mode contributions change with λ/λ_cr
lambda_norm = linspace(0, 1.5, 100);
% Smooth transition from Mode 1 to Mode 2
contribution_mode1 = cos(lambda_norm * pi/2);
contribution_mode2 = sin(lambda_norm * pi/2);

hold on;
% Filled areas
fill([lambda_norm, fliplr(lambda_norm)], [contribution_mode1, zeros(size(contribution_mode1))], ...
     'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([lambda_norm, fliplr(lambda_norm)], [zeros(size(contribution_mode2)), contribution_mode2], ...
     'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% Lines
plot(lambda_norm, contribution_mode1, 'b-', 'LineWidth', 2);
plot(lambda_norm, contribution_mode2, 'r-', 'LineWidth', 2);

% Mark flutter point (λ = λ_cr)
xline(1, 'k--', 'LineWidth', 2);
plot(1, 0.5, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'k');
text(1.02, 0.52, 'FLUTTER', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');

xlabel('λ / λ_{cr}', 'FontSize', 10);
ylabel('Modal Contribution', 'FontSize', 10);
title('Mode Coalescence at Flutter (CCCC)', 'FontSize', 11);
legend('Mode 1 (Bending)', 'Mode 2 (Torsion)', 'Location', 'east');
grid on;
xlim([0, 1.5]);
ylim([0, 1]);

% ========================================================================
% ROW 2, COL 2: MODE SHAPE EVOLUTION WITH λ (Centerline)
% ========================================================================
subplot(2,3,5);
% Show centerline deflection as λ increases
x_center = linspace(0, a, n_plot);
y_center = b/2;

% Different λ stages (normalized by λ_cr)
lambda_stages = [0, 0.3, 0.6, 0.9, 1.0, 1.2];
colors = {'b', 'c', 'g', 'm', 'r', 'k'};
line_styles = {'-', '-', '-', '-', '-', '--'};

hold on;
for s = 1:length(lambda_stages)
    if lambda_stages(s) <= 1.0
        weight1 = cos(lambda_stages(s) * pi/2);
        weight2 = sin(lambda_stages(s) * pi/2);
    else
        weight1 = 0;
        weight2 = 1;
    end
    
    Z_centerline = zeros(1, n_plot);
    for i = 1:n_plot
        Z_centerline(i) = weight1 * Z_mode1(round(n_plot/2), i) + weight2 * Z_mode2(round(n_plot/2), i);
    end
    Z_centerline = Z_centerline / max(abs(Z_centerline));
    
    if lambda_stages(s) == 1.0
        plot(x_center, Z_centerline, line_styles{s}, 'Color', colors{s}, 'LineWidth', 3);
    else
        plot(x_center, Z_centerline, line_styles{s}, 'Color', colors{s}, 'LineWidth', 1.5);
    end
end

% Mark clamped ends (CCCC - both ends clamped)
xline(0, 'r-', 'LineWidth', 2);
xline(a, 'r-', 'LineWidth', 2);
text(0.02, 0.8, 'CLAMPED', 'FontSize', 8, 'Color', 'r');
text(a-0.08, 0.8, 'CLAMPED', 'FontSize', 8, 'Color', 'r');

xlabel('Length x (m)', 'FontSize', 10);
ylabel('Normalized Displacement', 'FontSize', 10);
title('Mode Shape Evolution with λ (Centerline)', 'FontSize', 11);
legend({'λ=0', 'λ=0.3λ_{cr}', 'λ=0.6λ_{cr}', 'λ=0.9λ_{cr}', 'λ=λ_{cr} (FLUTTER)', 'λ=1.2λ_{cr}'}, ...
       'Location', 'southeast', 'FontSize', 7);
grid on;
xlim([0, a]);
ylim([-1.2, 1.2]);

% ========================================================================
% ROW 2, COL 3: NODAL PATTERN AT FLUTTER
% ========================================================================
subplot(2,3,6);
% Create contour plot showing nodal lines at flutter
Z_nodal = Z_coalesced;
levels = linspace(-0.8, 0.8, 9);
contourf(X_plot, Y_plot, Z_nodal, levels, 'EdgeColor', 'none');
colormap(jet);
colorbar;
hold on;

% Plot zero contour (nodal line)
[C, h] = contour(X_plot, Y_plot, Z_nodal, [0 0], 'k-', 'LineWidth', 2.5);
clabel(C, h, 'FontSize', 10, 'Color', 'k');

% Mark all clamped edges (CCCC)
plot([0 a a 0 0], [0 0 b b 0], 'r-', 'LineWidth', 2);
text(a/2, -0.015, 'CLAMPED', 'FontSize', 8, 'Color', 'r', 'HorizontalAlignment', 'center');
text(a/2, b + 0.01, 'CLAMPED', 'FontSize', 8, 'Color', 'r', 'HorizontalAlignment', 'center');
text(-0.025, b/2, 'CLAMPED', 'Rotation', 90, 'FontSize', 8, 'Color', 'r');
text(a + 0.01, b/2, 'CLAMPED', 'Rotation', -90, 'FontSize', 8, 'Color', 'r');

xlabel('Length x (m)', 'FontSize', 10);
ylabel('Width y (m)', 'FontSize', 10);
if flutter_detected
    title(sprintf('Nodal Pattern at Flutter (λ_{cr}=%.0f)', lambda_cr), 'FontSize', 11);
else
    title('Nodal Pattern at Flutter (Theoretical)', 'FontSize', 11);
end
axis equal;
xlim([-0.03, a+0.03]);
ylim([-0.03, b+0.03]);

% ========================================================================
% MAIN TITLE
% ========================================================================
if flutter_detected
    sgtitle(sprintf('CCCC PANEL MODE SHAPES: Flutter at λ_{cr}=%.0f (M=%.2f, f=%.1f Hz)', ...
            lambda_cr, Mach_flutter, freq_flutter), 'FontSize', 14, 'FontWeight', 'bold');
else
    sgtitle('CCCC PANEL MODE SHAPES (Theoretical Flutter Shape)', 'FontSize', 14, 'FontWeight', 'bold');
end

% Save figures
saveas(gcf, 'mode_shapes_cccc.png');
saveas(gcf, 'mode_shapes_cccc.fig');
fprintf('\n✓ Mode shapes saved as: mode_shapes_cccc.png\n');
%----------------------------------------------------------------------------
figure('Position', [100, 100, 1400, 500]);

% Subplot 1: Real part of eigenvalues
subplot(1,2,1);
if exist('real_eigenvalues', 'var') && ~isempty(real_eigenvalues)
    plot(lambda_values, real_eigenvalues(1,:), 'b-', 'LineWidth', 2); hold on;
    plot(lambda_values, real_eigenvalues(2,:), 'r-', 'LineWidth', 2);
else
    % If real_eigenvalues doesn't exist, plot zeros
    plot(lambda_values, zeros(size(lambda_values)), 'b-', 'LineWidth', 2); hold on;
    plot(lambda_values, zeros(size(lambda_values)), 'r-', 'LineWidth', 2);
    fprintf('Warning: real_eigenvalues not available\n');
end

if flutter_detected
    plot(lambda_cr, -18, 'ko', 'MarkerSize', 10, 'LineWidth', 2);
    xline(lambda_cr, 'k--', 'LineWidth', 1.5);
end
yline(0, 'k--', 'LineWidth', 1);
xlabel('Non-dimensional Aerodynamic Pressure \lambda', 'FontSize', 12);
ylabel('Real Part of Eigenvalue', 'FontSize', 12);
title('Real Part of Eigenvalues vs \lambda', 'FontSize', 14);
legend('Mode 1', 'Mode 2', 'Location', 'best');
grid on;
xlim([lambda_min, lambda_max]);

% Subplot 2: Damping ratio
subplot(1,2,2);
plot(lambda_values, damping_ratio(1,:), 'b-', 'LineWidth', 2); hold on;
plot(lambda_values, damping_ratio(2,:), 'r-', 'LineWidth', 2);
if flutter_detected
    plot(lambda_cr, 0, 'ko', 'MarkerSize', 10, 'LineWidth', 2);
    xline(lambda_cr, 'k--', 'LineWidth', 1.5);
end
yline(0, 'k--', 'LineWidth', 1);
xlabel('Non-dimensional Aerodynamic Pressure \lambda', 'FontSize', 12);
ylabel('Damping Ratio', 'FontSize', 12);
title('Damping Ratio vs \lambda', 'FontSize', 14);
legend('Mode 1', 'Mode 2', 'Location', 'best');
grid on;
xlim([lambda_min, lambda_max]);

if flutter_detected
    sgtitle(sprintf('Flutter Analysis: λ_{cr} = %.1f, f = %.0f Hz (Mach %.2f)', ...
            lambda_cr, freq_flutter, Mach_flutter), 'FontSize', 14);
else
    sgtitle('Flutter Analysis: No Flutter Detected', 'FontSize', 14);
end


%% ========================================================================
% SECTION 8: GAUSS-LEGENDRE FUNCTION
% ========================================================================
function [x, w] = gauss_legendre(n, a, b)
    beta = 0.5 ./ sqrt(1 - (2*(1:n-1)).^(-2));
    T = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(T);
    x_leg = diag(D);
    w_leg = 2 * V(1,:)'.^2;
    [x_leg, idx] = sort(x_leg);
    w_leg = w_leg(idx);
    x = (b - a)/2 * x_leg + (a + b)/2;
    w = (b - a)/2 * w_leg;
end

%% ========================================================================
% SECTION 9: LAMINATE PROPERTIES CALCULATION
% ========================================================================
function [A, B, D, As, I0, I2] = calculate_laminate_properties(layers)
    if nargin < 1
        error('At least one input argument (layers) is required');
    end
    if ~iscell(layers)
        error('layers must be a cell array');
    end
    
    A = zeros(3);
    B = zeros(3);
    D = zeros(3);
    As = zeros(2);
    I0 = 0;
    I2 = 0;
    
    h_tot = sum(cellfun(@(x) x.thickness, layers));
    z_cur = -h_tot / 2;
    
    for i = 1:length(layers)
        L = layers{i};
        
        if ~isfield(L, 'type') || ~isfield(L, 'thickness') || ~isfield(L, 'rho')
            error('Layer %d is missing required fields', i);
        end
        
        z_next = z_cur + L.thickness;
        
        if strcmp(L.type, 'comp')
            required_fields = {'E1', 'E2', 'G12', 'nu12', 'angle'};
            for f = 1:length(required_fields)
                if ~isfield(L, required_fields{f})
                    error('Layer %d missing field: %s', i, required_fields{f});
                end
            end
            
            nu21 = L.nu12 * L.E2 / L.E1;
            denom = 1 - L.nu12 * nu21;
            Q11 = L.E1 / denom;
            Q22 = L.E2 / denom;
            Q12 = L.nu12 * L.E2 / denom;
            Q66 = L.G12;
            
            Q_local = [Q11, Q12, 0; Q12, Q22, 0; 0, 0, Q66];
            
            theta = L.angle * pi / 180;
            m = cos(theta);
            n = sin(theta);
            
            T = [m^2, n^2, 2*m*n; n^2, m^2, -2*m*n; -m*n, m*n, m^2-n^2];
            Q_trans = T \ Q_local * T;
            
            Qs_local = [L.G12, 0; 0, L.G12];
            T_shear = [m^2, n^2; n^2, m^2];
            Qs_trans = T_shear * Qs_local;
            
        elseif strcmp(L.type, 'piezo')
            if ~isfield(L, 'E') || ~isfield(L, 'nu')
                error('Layer %d missing field: E or nu', i);
            end
            
            E = L.E;
            nu = L.nu;
            denom = 1 - nu^2;
            Q_trans = (E / denom) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
            Qs_trans = (E / (2 * (1 + nu))) * eye(2);
        else
            error('Layer %d unknown type: %s', i, L.type);
        end
        
        dz = z_next - z_cur;
        dz2 = z_next^2 - z_cur^2;
        dz3 = z_next^3 - z_cur^3;
        
        A = A + Q_trans * dz;
        B = B + Q_trans * dz2 / 2;
        D = D + Q_trans * dz3 / 3;
        As = As + Qs_trans * dz;
        
        I0 = I0 + L.rho * dz;
        I2 = I2 + L.rho * dz3 / 3;
        
        z_cur = z_next;
    end
    
    As = As * 5/6;
    
    if norm(B) < 1e-10
        B = zeros(3);
    end
    
    A = (A + A') / 2;
    B = (B + B') / 2;
    D = (D + D') / 2;
    As = (As + As') / 2;
end
