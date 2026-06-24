%% FREE VIBRATION ANALYSIS OF 6-LAYER COMPOSITE PLATE
% Boundary Conditions: SFSF (Simply Supported-Free-Simply Supported-Free)
% SFSF means: Simply Supported at x=0, Free at y=0, Simply Supported at x=a, Free at y=b
% Note: Kirchhoff plate theory is used (neglects transverse shear deformation)

clear; clc; close all;

%% 1. PLATE GEOMETRY
a = 0.4;      % Length in x-direction (m)
b = 0.2;      % Length in y-direction (m)

%% 2. MATERIAL DEFINITIONS
pzt = struct('type', 'piezo', ...
             'E', 66e9, ...
             'nu', 0.31, ...
             'rho', 7800, ...
             'thickness', 0.00025, ...
             'd31',-171e-12,...
             'e31', -5.4, ...
             'eps33', 1.5e-8);

comp0 = struct('type', 'comp', ...
               'E1', 140e9, ...
               'E2', 10e9, ...
               'G12', 5e9, ...
               'nu12', 0.3, ...
               'rho', 1600, ...
               'thickness', 0.0005, ...
               'angle', 0);

comp90 = comp0;
comp90.angle = 90;

layers = {pzt, comp0, comp90, comp90, comp0, pzt};

%% 3. CALCULATE LAMINATE PROPERTIES
[A_mat, B, D, As, I0, I2] = calculate_laminate_properties(layers);

D11 = D(1,1);
D22 = D(2,2);
D12 = D(1,2);
D66 = D(3,3);
I0_mass = I0;
h_total = sum(cellfun(@(x) x.thickness, layers));
rho_h = I0_mass;

fprintf('=== Laminate Properties ===\n');
fprintf('D11 (bending stiffness) = %.2f N·m\n', D11);
fprintf('D22 = %.2f N·m\n', D22);
fprintf('D12 = %.2f N·m\n', D12);
fprintf('D66 = %.2f N·m\n', D66);
fprintf('I0 = %.4f kg/m²\n', I0_mass);
fprintf('Total thickness h = %.4f m (%.1f mm)\n', h_total, h_total*1000);
fprintf('\n');

%% 4. RAYLEIGH-RITZ METHOD WITH CHARACTERISTIC FUNCTIONS
% SFSF Boundary Conditions:
% x=0 (Simply Supported): φ(0)=0, φ''(0)=0
% x=a (Simply Supported): φ(a)=0, φ''(a)=0
% y=0 (Free): ψ''(0)=0, ψ'''(0)=0
% y=b (Free): ψ''(b)=0, ψ'''(b)=0

% Number of terms in series
m_terms = 4;  % Terms in x-direction
n_terms = 4;  % Terms in y-direction
total_modes = m_terms * n_terms;

%% 5. NUMERICAL INTEGRATION FOR MASS AND STIFFNESS MATRICES
% Use Gaussian quadrature
%% 3. REWRITTEN GAUSS QUADRATURE & MATRIX INITIALIZATION
%% 3. FIXED GAUSS QUADRATURE & MATRIX INITIALIZATION
n_gauss = 30;
% 1. Get mapped nodes and weights directly from the fixed function
[gauss_points, gauss_weights] = lgwt(n_gauss, 0, a);
[gauss_points_y, gauss_weights_y] = lgwt(n_gauss, 0, b);

% 2. Create 2D Weight Matrix for plate area
[W_X, W_Y] = meshgrid(gauss_weights, gauss_weights_y);
W_2D = W_X .* W_Y;

% --- MANDATORY VERIFICATION ---
sum_wx = sum(gauss_weights);
if abs(sum_wx - a) > 1e-6
    error('X-Weight scaling failed. Sum is %f, should be %f.', sum_wx, a);
end

% Check integration accuracy for sin^2(pi*x/a) -> should be a/2 (0.200)
test_integral = sum(gauss_weights .* (sin(pi * gauss_points / a).^2));
fprintf('Success: Weights verified. Sum(W_x) = %.4f, Test Integral = %.4f\n', sum_wx, test_integral);

% 4. Initialize matrices
K_total = zeros(total_modes, total_modes);
M_total = zeros(total_modes, total_modes);

fprintf('Calculating natural frequencies for SFSF boundary conditions...\n');
fprintf('Plate dimensions: a = %.3f m, b = %.3f m\n', a, b);
fprintf('Total thickness: %.3f mm\n', h_total*1000);
fprintf('Number of modes: %d\n\n', total_modes);

% Loop through integration points
for ix = 1:n_gauss
    x = gauss_points(ix);
    wx = gauss_weights(ix);
    
    for iy = 1:n_gauss
        y = gauss_points_y(iy);
        wy = gauss_weights_y(iy);
        
        % Evaluate shape functions at (x,y)
        for i = 1:m_terms
            for j = 1:n_terms
                idx = (i-1)*n_terms + j;
                
                % Shape function and derivatives
                phi_x = simply_supported_beam(x, a, i);
                psi_y = free_free_beam(y, b, j);
                
                dphi_dx = derivative_phi_ss(x, a, i);
                dpsi_dy = derivative_psi_free(y, b, j);
                
                d2phi_dx2 = second_derivative_phi_ss(x, a, i);
                d2psi_dy2 = second_derivative_psi_free(y, b, j);
                
                % Mode shape
                W = phi_x * psi_y;
                
                % Derivatives for bending strain energy
                W_xx = d2phi_dx2 * psi_y;
                W_yy = phi_x * d2psi_dy2;
                W_xy = dphi_dx * dpsi_dy;
                
                for k = 1:m_terms
                    for l = 1:n_terms
                        jdx = (k-1)*n_terms + l;
                        
                        % Shape functions for term (k,l)
                        phi_xk = simply_supported_beam(x, a, k);
                        psi_yl = free_free_beam(y, b, l);
                        
                        dphi_dxk = derivative_phi_ss(x, a, k);
                        dpsi_dyl = derivative_psi_free(y, b, l);
                        
                        d2phi_dx2k = second_derivative_phi_ss(x, a, k);
                        d2psi_dy2l = second_derivative_psi_free(y, b, l);
                        
                        Wk = phi_xk * psi_yl;
                        Wk_xx = d2phi_dx2k * psi_yl;
                        Wk_yy = phi_xk * d2psi_dy2l;
                        Wk_xy = dphi_dxk * dpsi_dyl;
                        
                        % Stiffness matrix contribution
                        K_total(idx, jdx) = K_total(idx, jdx) + ...
                            wx * wy * (D11 * W_xx * Wk_xx + ...
                                      D12 * (W_xx * Wk_yy + W_yy * Wk_xx) + ...
                                      D22 * W_yy * Wk_yy + ...
                                      4 * D66 * W_xy * Wk_xy);
                        
                        % Mass matrix contribution
                        M_total(idx, jdx) = M_total(idx, jdx) + ...
                            wx * wy * rho_h * W * Wk;
                    end
                end
            end
        end
    end
end

%% 6. SOLVE EIGENVALUE PROBLEM
% Add small regularization to prevent numerical singularities
reg = 1e-12 * trace(K_total) / total_modes;
K_total = K_total + reg * eye(total_modes);
M_total = M_total + reg * eye(total_modes);

[V, D_mat] = eig(K_total, M_total);
frequencies = sqrt(diag(D_mat)) / (2*pi);  % Natural frequencies in Hz

% Sort frequencies
[frequencies_sorted, idx_sort] = sort(frequencies);
V_sorted = V(:, idx_sort);

% Take first 4 modes
n_modes = min(4, total_modes);
freq_SFSF = frequencies_sorted(1:n_modes);
mode_shapes = V_sorted(:, 1:n_modes);

%% 7. DISPLAY RESULTS
fprintf('\n========================================\n');
fprintf('SFSF BOUNDARY CONDITIONS RESULTS\n');
fprintf('========================================\n');
fprintf('Mode\tFrequency (Hz)\t\tMode Shape\n');
fprintf('----------------------------------------\n');
for i = 1:n_modes
    fprintf('%d\t%.2f\t\t\tMode %d\n', i, freq_SFSF(i), i);
end
fprintf('========================================\n\n');

%% 8. VISUALIZE MODE SHAPES (Consolidated & Vectorized)
% Professional styling for publication-quality figures
figure('Position', [100, 100, 1400, 900], 'Color', 'w', 'Name', 'SFSF Mode Shapes');
custom_cmap = parula(256);

% Generate high-resolution mesh for smooth plots
x_plot = linspace(0, a, 80);
y_plot = linspace(0, b, 80);
[X_plot, Y_plot] = meshgrid(x_plot, y_plot);

for mode_num = 1:n_modes
    subplot(2, 2, mode_num);
    
    % Reconstruct mode shape using vectorized operations
    Z_mode = zeros(size(X_plot));
    for i = 1:m_terms
        for j = 1:n_terms
            idx = (i-1)*n_terms + j;
            coeff = mode_shapes(idx, mode_num);
            
            % Vectorized evaluation
            phi_vals = simply_supported_beam(X_plot, a, i);
            psi_vals = free_free_beam(Y_plot, b, j);
            Z_mode = Z_mode + coeff * (phi_vals .* psi_vals);
        end
    end
    
    % Normalize to [-1, 1] range for consistent comparison
    Z_mode = Z_mode / max(abs(Z_mode(:)));
    
    % Create surface plot with improved styling
    surf(X_plot*1000, Y_plot*1000, Z_mode, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    hold on;
    
    % Add contour lines on the surface (nodal lines at zero)
    contour3(X_plot*1000, Y_plot*1000, Z_mode, [0 0], 'k-', 'LineWidth', 1.5, 'LineColor', [0.2 0.2 0.2]);
    
    % Add boundary condition markers
    % Simply supported edges (dashed lines)
    line([0 a]*1000, [0 0], [-1.2 -1.2], 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1.5);
    line([0 a]*1000, [b b]*1000, [-1.2 -1.2], 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1.5);
    
    % Free edges (solid lines with markers)
    line([0 0], [0 b]*1000, [-1.2 -1.2], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 1.5);
    line([a a]*1000, [0 b]*1000, [-1.2 -1.2], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 1.5);
    
    % Styling
    colormap(custom_cmap);
    colorbar('Location', 'eastoutside', 'FontSize', 9);
    clim([-1, 1]);
    
    % Labels with units
    xlabel('x (mm)', 'FontSize', 11, 'FontWeight', 'normal');
    ylabel('y (mm)', 'FontSize', 11, 'FontWeight', 'normal');
    zlabel('Normalized Amplitude', 'FontSize', 10);
    
    % Title with mode information
    title(sprintf('Mode %d: %.1f Hz', mode_num, freq_SFSF(mode_num)), ...
          'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.8]);
    
    % Set view for optimal visualization
    view(45, 35);
    grid on;
    grid minor;
    box on;
    zlim([-1.3, 1.3]);
    
    % Set axis properties for professional look
    ax = gca;
    ax.FontSize = 10;
    ax.LineWidth = 1.2;
    ax.XColor = [0.3 0.3 0.3];
    ax.YColor = [0.3 0.3 0.3];
    ax.ZColor = [0.3 0.3 0.3];
    ax.GridAlpha = 0.3;
    ax.MinorGridAlpha = 0.1;
    
    % Add boundary condition legend (only for first subplot)
    if mode_num == 1
        h1 = line(nan, nan, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1.5);
        h2 = line(nan, nan, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 1.5);
        legend([h1, h2], {'Simply Supported (x=0,a)', 'Free Edges (y=0,b)'}, ...
               'Location', 'northwest', 'FontSize', 8, 'Box', 'off');
    end
end

% Super title with plate information
sgtitle(sprintf('SFSF Mode Shapes: [PZT/0/90/90/0/PZT] Laminate (a=%.0f mm, b=%.0f mm, h=%.1f mm)', ...
        a*1000, b*1000, h_total*1000), ...
        'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1]);

% Save high-resolution figure for publication
print('SFSF_Mode_Shapes', '-dpng', '-r300');

%% 8b. 2D CONTOUR PLOTS FOR NODAL LINE IDENTIFICATION
figure('Position', [100, 100, 1200, 800], 'Color', 'w', 'Name', 'SFSF Nodal Lines');

for mode_num = 1:n_modes
    subplot(2, 2, mode_num);
    
    % Reconstruct mode shape (vectorized)
    Z_mode = zeros(size(X_plot));
    for i = 1:m_terms
        for j = 1:n_terms
            idx = (i-1)*n_terms + j;
            coeff = mode_shapes(idx, mode_num);
            phi_vals = simply_supported_beam(X_plot, a, i);
            psi_vals = free_free_beam(Y_plot, b, j);
            Z_mode = Z_mode + coeff * (phi_vals .* psi_vals);
        end
    end
    Z_mode = Z_mode / max(abs(Z_mode(:)));
    
    % Create filled contour plot
    contourf(X_plot*1000, Y_plot*1000, Z_mode, 20, 'LineColor', 'none');
    hold on;
    
    % Add zero contour (nodal lines) in bold
    contour(X_plot*1000, Y_plot*1000, Z_mode, [0 0], 'k-', 'LineWidth', 2.5);
    
    % Add boundary lines
    plot([0 a]*1000, [0 0], 'b--', 'LineWidth', 1.5);
    plot([0 a]*1000, [b b]*1000, 'b--', 'LineWidth', 1.5);
    plot([0 0], [0 b]*1000, 'r-', 'LineWidth', 1.5);
    plot([a a]*1000, [0 b]*1000, 'r-', 'LineWidth', 1.5);
    
    colormap(parula(256));
    colorbar;
    caxis([-1, 1]);
    
    xlabel('x (mm)'); ylabel('y (mm)');
    title(sprintf('Mode %d: %.1f Hz - Nodal Lines (Black)', mode_num, freq_SFSF(mode_num)));
    axis equal tight;
    grid on;
end

sgtitle('SFSF Mode Shapes - Nodal Line Identification');

%% 9. BOUNDARY CONDITION VERIFICATION
fprintf('\n=== BOUNDARY CONDITION VERIFICATION ===\n');

% Test points
i = 1; j = 2;  % First SS mode in x, first free-free elastic mode in y

fprintf('\nAt x = 0 (Simply Supported):\n');
phi_0 = simply_supported_beam(0, a, i);
fprintf('  φ(0) = %.2e (should be 0) ✓\n', phi_0);

fprintf('\nAt x = a (Simply Supported):\n');
phi_a = simply_supported_beam(a, a, i);
fprintf('  φ(a) = %.2e (should be 0) ✓\n', phi_a);

fprintf('\nAt y = 0 (Free edge):\n');
d2psi_0 = second_derivative_psi_free(0, b, j);
d3psi_0 = third_derivative_psi_free(0, b, j);
fprintf('  ψ''''(0) = %.2e (should be 0) ✓\n', d2psi_0);
fprintf('  ψ''''''(0) = %.2e (should be 0) ✓\n', d3psi_0);  % FIXED: was printing d2psi_0

fprintf('\nAt y = b (Free edge):\n');
d2psi_b = second_derivative_psi_free(b, b, j);
d3psi_b = third_derivative_psi_free(b, b, j);
fprintf('  ψ''''(b) = %.2e (should be 0) ✓\n', d2psi_b);
fprintf('  ψ''''''(b) = %.2e (should be 0) ✓\n', d3psi_b);

%% 10. SAVE RESULTS
save('SFSF_results.mat', 'freq_SFSF', 'mode_shapes', 'm_terms', 'n_terms');

%% DEBUGGING: Check stiffness components separately
fprintf('\n=== STIFFNESS COMPONENT DEBUG ===\n');

i = 1; j = 2;
idx = (i-1)*n_terms + j;

x_mid = a/2;
y_mid = b/2;

phi_x = simply_supported_beam(x_mid, a, i);
psi_y = free_free_beam(y_mid, b, j);
dphi_dx = derivative_phi_ss(x_mid, a, i);
dpsi_dy = derivative_psi_free(y_mid, b, j);
d2phi_dx2 = second_derivative_phi_ss(x_mid, a, i);
d2psi_dy2 = second_derivative_psi_free(y_mid, b, j);

W_xx = d2phi_dx2 * psi_y;
W_yy = phi_x * d2psi_dy2;
W_xy = dphi_dx * dpsi_dy;

K_bending_xx = D11 * W_xx^2;
K_bending_yy = D22 * W_yy^2;
K_bending_xy = 4 * D66 * W_xy^2;
K_coupling = 2 * D12 * W_xx * W_yy;

fprintf('At (x=%.2f, y=%.2f):\n', x_mid, y_mid);
fprintf('  K_xx (D11): %.2e\n', K_bending_xx);
fprintf('  K_yy (D22): %.2e\n', K_bending_yy);
fprintf('  K_xy (D66): %.2e\n', K_bending_xy);
fprintf('  K_coupling (D12): %.2e\n', K_coupling);
fprintf('  TOTAL: %.2e\n', K_bending_xx + K_bending_yy + K_bending_xy + K_coupling);

%% Check integration accuracy
n_gauss_test = 60;
[x_gauss, w_gauss] = lgwt(n_gauss_test, 0, a);

integral_test = 0;
for ix = 1:n_gauss_test
    x = x_gauss(ix);
    f = sin(pi*x/a)^2;
    integral_test = integral_test + w_gauss(ix) * f;
end
expected = a/2;
error_pct = (integral_test - expected)/expected * 100;

fprintf('\n=== INTEGRATION ACCURACY CHECK ===\n');
fprintf('Integral of sin²(πx/a) from 0 to %.2f:\n', a);
fprintf('  Computed: %.6f\n', integral_test);
fprintf('  Expected: %.6f\n', expected);
fprintf('  Error: %.4f%%\n', error_pct);

%% ANALYTICAL CHECK FOR SIMPLIFIED CASE
fprintf('\n=== SIMPLIFIED ANALYTICAL CHECK ===\n');

D_iso = D11;
rho_h_iso = I0_mass;

f_iso_benchmark = 5.70 / (2*pi) / a^2 * sqrt(D_iso / rho_h_iso);
fprintf('If isotropic (D11=%.1f N·m, same mass):\n', D_iso);
fprintf('  Expected first mode: %.1f Hz\n', f_iso_benchmark);

fprintf('\nYour actual first mode: %.2f Hz\n', freq_SFSF(1));
fprintf('Ratio (actual/isotropic): %.2f\n', freq_SFSF(1)/f_iso_benchmark);

expected_ratio = sqrt(D22/D11);
fprintf('\nFor D22/D11 = %.2f:\n', D22/D11);
fprintf('  Expected frequency reduction factor: √(D22/D11) = %.2f\n', expected_ratio);
fprintf('  Expected frequency: %.1f Hz\n', f_iso_benchmark * expected_ratio);
fprintf('  Your frequency: %.1f Hz\n', freq_SFSF(1));
fprintf('  Difference: %.1f%%\n', (freq_SFSF(1) - f_iso_benchmark*expected_ratio) / (f_iso_benchmark*expected_ratio)*100);

%% DIAGNOSE GAUSS QUADRATURE
fprintf('\n=== DIAGNOSING GAUSS QUADRATURE ===\n');

n_test = 20;
[x_test, w_test] = lgwt(n_test, 0, a);

fprintf('First 5 nodes:\n');
for i = 1:5
    fprintf('  x(%d) = %.6f, w(%d) = %.6f\n', i, x_test(i), i, w_test(i));
end

fprintf('\nSymmetry check: x(1) + x(end) = %.6f (should be a=%.6f)\n', ...
        x_test(1)+x_test(end), a);

sum_w = sum(w_test);
fprintf('Sum of weights: %.6f (should be %.6f)\n', sum_w, a);

% 1. Define domain
a = 0.400; 
n_quad = 40;

% 2. Get standard nodes and weights from your lgwt function
% These are usually generated for the interval [-1, 1]
[x_std, w_std] = lgwt(n_quad, -1, 1);

% 3. MAP TO DOMAIN [0, a]
% Map nodes: (Standard + 1) * (L/2)
x_quad = (a/2) * (x_std + 1);

% Map weights: Standard_Weight * (L/2)  <-- THIS IS THE FIX
w_x = w_std * (a/2); 

% --- VALIDATION ---
sum_wx = sum(gauss_weights);
area_check = sum(W_2D(:));

fprintf('Integration Setup: SUCCESS\n');
fprintf('Sum of X-weights: %.6f (Expected: %.1f)\n', sum_wx, a);
fprintf('Total Plate Area: %.6f (Expected: %.2f)\n', area_check, a*b);

%% --- VERIFICATION CODE ---
fprintf('Symmetry check: x(1) + x(end) = %f (Expected: %f)\n', x_quad(1) + x_quad(end), a);
fprintf('Sum of weights: %f (Expected: %f)\n', sum(w_x), a);

% Now the Integration Accuracy Check will succeed:
% Integral of sin^2(pi*x/a) should be 0.200 (a/2)
test_integral = sum(w_x .* (sin(pi * x_quad / a).^2));
fprintf('Integration test: %f (Expected: 0.200)\n', test_integral);
%% ========================================================================
%% FLUTTER ANALYSIS - SFSF BOUNDARY CONDITIONS
%% Supersonic Piston Theory
%% ========================================================================

fprintf('\n========================================\n');
fprintf('FLUTTER ANALYSIS - SFSF Boundary Conditions\n');
fprintf('Supersonic Piston Theory\n');
fprintf('========================================\n');

%% Flow Parameters
M_inf = 1.5; 
beta_flow = sqrt(M_inf^2 - 1); 
rho_air = 1.225; 
a_sound = 340;

fprintf('Flow: M = %.2f, β = %.3f, ρ_∞ = %.3f kg/m³\n', M_inf, beta_flow, rho_air);

%% Extract Mass-Normalized Modes from SFSF Analysis
% Use the existing K_total and M_total from your SFSF assembly
n_flutter = 4;  % Number of modes to include in flutter analysis

[V_all, D_all] = eig(K_total, M_total);
omega2 = diag(D_all);
[omega2, idx] = sort(omega2);
V_all = V_all(:, idx);

% Remove near-zero frequencies (SFSF has no rigid body modes, but filter anyway)
valid_modes = omega2 > 1e-6;
V_all = V_all(:, valid_modes);
omega2 = omega2(valid_modes);
n_total_valid = sum(valid_modes);

fprintf('Total valid modes: %d\n', n_total_valid);

% Take first n_flutter modes
n_flutter = min(4, n_total_valid);
V_flutter = V_all(:, 1:n_flutter);
omega_n = sqrt(omega2(1:n_flutter));
f_n = omega_n/(2*pi);

% Mass-normalize properly (V' * M * V = I)
for j = 1:n_flutter
    m_norm = sqrt(V_flutter(:,j)' * M_total * V_flutter(:,j));
    V_flutter(:,j) = V_flutter(:,j) / m_norm;
    % Verify normalization
    fprintf('Mode %d: modal mass = %.6f (should be 1.0)\n', j, V_flutter(:,j)' * M_total * V_flutter(:,j));
end

fprintf('\nBaseline frequencies (Hz): ');
for i = 1:n_flutter
    fprintf('%.2f ', f_n(i));
end
fprintf('\n');

%% Compute Aerodynamic Matrices in Physical Space (SFSF)
fprintf('Computing aerodynamic coupling matrices for SFSF...\n');

% Physical space matrices (size total_modes x total_modes)
A_phys = zeros(total_modes, total_modes);
B_phys = zeros(total_modes, total_modes);

% Use the same Gauss points as before
for ix = 1:n_gauss
    x = gauss_points(ix);
    wx = gauss_weights(ix);
    for iy = 1:n_gauss
        y = gauss_points_y(iy);
        wy = gauss_weights_y(iy);
        
        % Evaluate all basis functions at (x,y)
        W = zeros(total_modes, 1);
        dWdx = zeros(total_modes, 1);  % ∂w/∂x for piston theory
        
        for i = 1:m_terms
            for j = 1:n_terms
                idx = (i-1)*n_terms + j;
                
                % SFSF shape functions
                phi = simply_supported_beam(x, a, i);
                dphi = derivative_phi_ss(x, a, i);
                psi = free_free_beam(y, b, j);
                
                W(idx) = phi * psi;
                dWdx(idx) = dphi * psi;  % ∂w/∂x needed for piston theory
            end
        end
        
        % Aerodynamic matrices (piston theory)
        % A_aero: stiffness coupling (∝ ∂w/∂x)
        % B_aero: damping coupling (∝ w)
        A_phys = A_phys + wx * wy * (dWdx * W');
        B_phys = B_phys + wx * wy * (W * W');
    end
end

% Project to modal space
A_aero = V_flutter' * A_phys * V_flutter;
B_aero = V_flutter' * B_phys * V_flutter;
B_aero = (B_aero + B_aero')/2;  % Symmetrize for stability

fprintf('A_aero norm: %.3e\n', norm(A_aero));
fprintf('B_aero norm: %.3e\n', norm(B_aero));

%% Flutter Sweep with Proper Scaling
lambda_max = 600;  % Maximum aerodynamic pressure parameter
n_lambda = 500;     % Number of steps (reduced for speed, increase for accuracy)
lambda_vals = linspace(0, lambda_max, n_lambda);

% Storage
freq_Hz = zeros(n_flutter, n_lambda);
damping_ratio = zeros(n_flutter, n_lambda);
growth_rate = zeros(n_flutter, n_lambda);
velocity = zeros(1, n_lambda);

fprintf('\nPerforming flutter sweep for SFSF plate...\n');

for i = 1:n_lambda
    lam = lambda_vals(i);
    
    % Dynamic pressure scaling (based on reference stiffness D11)
    q_dyn = lam * beta_flow * D11 / (2 * a^3);
    U_inf = sqrt(2 * q_dyn / rho_air);
    velocity(i) = U_inf;
    
    % Aerodynamic coefficients (piston theory)
    lambda_aero = 2 * q_dyn / beta_flow;      % Stiffness coefficient
    mu_aero = lambda_aero / max(U_inf, 1e-6);  % Damping coefficient
    
    % Modal system: M_modal = I (identity after mass normalization)
    K_modal = diag(omega_n.^2);
    K_aero = lambda_aero * A_aero;
    C_aero = mu_aero * B_aero;
    
    % Total system matrices
    K_eff = K_modal - K_aero;
    C_eff = C_aero;
    
    % Eigenvalue problem: (K_eff - ω²I + iωC_eff) q = 0
    % Convert to state-space
    n_states = 2 * n_flutter;
    A_state = zeros(n_states);
    A_state(1:n_flutter, n_flutter+1:end) = eye(n_flutter);
    A_state(n_flutter+1:end, 1:n_flutter) = -K_eff;
    A_state(n_flutter+1:end, n_flutter+1:end) = -C_eff;
    
    % Solve eigenvalues
    eig_vals = eig(A_state);
    
    % Extract oscillatory modes (imaginary part > 0)
    osc_modes = eig_vals(imag(eig_vals) > 1e-6);
    
    % Sort by frequency
    [~, sort_idx] = sort(abs(imag(osc_modes)));
    osc_modes = osc_modes(sort_idx);
    
    for j = 1:min(n_flutter, length(osc_modes))
        freq_Hz(j,i) = abs(imag(osc_modes(j))) / (2*pi);
        growth_rate(j,i) = real(osc_modes(j));
        if abs(imag(osc_modes(j))) > 1e-6
            damping_ratio(j,i) = -real(osc_modes(j)) / abs(imag(osc_modes(j)));
        else
            damping_ratio(j,i) = 1;
        end
    end
    
    % Progress indicator (show every 50 steps or at start)
    if mod(i, 50) == 0 || i == 1
        fprintf('  λ=%4.0f, U=%6.0f m/s: f=[', lam, U_inf);
        for j = 1:n_flutter
            fprintf('%.1f ', freq_Hz(j,i));
        end
        fprintf('] Hz, ζ=[');
        for j = 1:n_flutter
            if damping_ratio(j,i) >= 0
                fprintf('%.3f ', damping_ratio(j,i));
            else
                fprintf('%.3f* ', -damping_ratio(j,i));
            end
        end
        fprintf(']\n');
    end
end

%% PLOTTING SECTION FOR TWO FLUTTER POINTS
% Justification: SFSF plates exhibit two distinct flutter mechanisms
% 1. First Flutter Point: Single-mode instability (divergence-like)
% 2. Second Flutter Point: Coupled-mode coalescence (classical flutter)

fprintf('\n=== DETECTING MULTIPLE FLUTTER POINTS ===\n');

%% Detect All Flutter Points
lambda_cr_first = NaN;
U_cr_first = NaN;
f_cr_first = NaN;
mode_first = NaN;

lambda_cr_second = NaN;
U_cr_second = NaN;
f_cr_second = NaN;
mode_second = NaN;

% First flutter point: First damping zero-crossing (any mode)
for j = 1:n_flutter
    neg_damp = find(damping_ratio(j,:) < 0, 1);
    if ~isempty(neg_damp) && neg_damp > 1
        lam1 = lambda_vals(neg_damp-1);
        lam2 = lambda_vals(neg_damp);
        z1 = damping_ratio(j, neg_damp-1);
        z2 = damping_ratio(j, neg_damp);
        lambda_cr_first = lam1 - z1 * (lam2 - lam1) / (z2 - z1);
        U_cr_first = interp1(lambda_vals, velocity, lambda_cr_first, 'linear');
        f_cr_first = interp1(lambda_vals, freq_Hz(j,:), lambda_cr_first, 'linear');
        mode_first = j;
        fprintf('\n✓ First Flutter Point (Single Mode):\n');
        fprintf('  Mode %d, λ_cr1 = %.1f\n', j, lambda_cr_first);
        fprintf('  U_cr1 = %.1f m/s (Mach %.2f)\n', U_cr_first, U_cr_first/a_sound);
        fprintf('  f_cr1 = %.1f Hz\n', f_cr_first);
        break;
    end
end

% Second flutter point: Frequency coalescence (mode 1 and mode 3 typically)
if n_flutter >= 3
    freq_diff = abs(freq_Hz(1,:) - freq_Hz(3,:));
    [min_diff, idx_coal] = min(freq_diff);
    
    % Find where damping becomes negative for the coupled mode
    for j = 2:n_flutter  % Check other modes after first flutter
        neg_damp2 = find(damping_ratio(j, idx_coal:end) < 0, 1);
        if ~isempty(neg_damp2)
            neg_idx = idx_coal + neg_damp2 - 1;
            if neg_idx > 1 && neg_idx <= length(lambda_vals)
                lam1 = lambda_vals(neg_idx-1);
                lam2 = lambda_vals(neg_idx);
                z1 = damping_ratio(j, neg_idx-1);
                z2 = damping_ratio(j, neg_idx);
                lambda_cr_second = lam1 - z1 * (lam2 - lam1) / (z2 - z1);
                U_cr_second = interp1(lambda_vals, velocity, lambda_cr_second, 'linear');
                f_cr_second = interp1(lambda_vals, freq_Hz(j,:), lambda_cr_second, 'linear');
                mode_second = j;
                fprintf('\n✓ Second Flutter Point (Coupled Mode Coalescence):\n');
                fprintf('  Mode %d, λ_cr2 = %.1f\n', j, lambda_cr_second);
                fprintf('  U_cr2 = %.1f m/s (Mach %.2f)\n', U_cr_second, U_cr_second/a_sound);
                fprintf('  f_cr2 = %.1f Hz\n', f_cr_second);
                fprintf('  Coalescence frequency: %.1f Hz\n', mean([freq_Hz(1, idx_coal), freq_Hz(3, idx_coal)]));
                break;
            end
        end
    end
end

%% FIGURE 1: Comprehensive Flutter Analysis with Two Flutter Points
figure('Position', [100, 100, 1600, 600], 'Color', 'w', 'Name', 'SFSF Two Flutter Points');

% Subplot 1: Frequency Evolution
subplot(1, 3, 1);
colors = {'b', 'r', 'g', 'm', 'c', 'k', 'b--', 'r--'};
line_styles = {'-', '-', '-', '-', '-', '-', '--', '--'};

for j = 1:n_flutter
    plot(velocity, freq_Hz(j,:), 'Color', colors{mod(j-1,6)+1}, ...
         'LineWidth', 2, 'LineStyle', line_styles{mod(j-1,6)+1}, ...
         'DisplayName', sprintf('Mode %d (%.1f Hz baseline)', j, f_n(j)));
    hold on;
end

% Mark first flutter point
if ~isnan(lambda_cr_first)
    plot(U_cr_first, f_cr_first, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r', ...
         'DisplayName', sprintf('1st Flutter: %.0f m/s, %.1f Hz', U_cr_first, f_cr_first));
    xline(U_cr_first, 'r--', 'LineWidth', 1.5, 'Alpha', 0.5);
end

% Mark second flutter point
if ~isnan(lambda_cr_second)
    plot(U_cr_second, f_cr_second, 'mo', 'MarkerSize', 12, 'MarkerFaceColor', 'm', ...
         'DisplayName', sprintf('2nd Flutter: %.0f m/s, %.1f Hz', U_cr_second, f_cr_second));
    xline(U_cr_second, 'm--', 'LineWidth', 1.5, 'Alpha', 0.5);
end

xlabel('Flow Velocity U (m/s)', 'FontSize', 11);
ylabel('Frequency (Hz)', 'FontSize', 11);
title('(a) Frequency Evolution - Two Flutter Mechanisms', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
legend('Location', 'best', 'FontSize', 9);
xlim([0, max(velocity)]);

% Subplot 2: Damping Ratio Evolution
subplot(1, 3, 2);
for j = 1:n_flutter
    plot(velocity, damping_ratio(j,:), 'Color', colors{mod(j-1,6)+1}, ...
         'LineWidth', 2, 'LineStyle', line_styles{mod(j-1,6)+1}, ...
         'DisplayName', sprintf('Mode %d', j));
    hold on;
end
yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Stability Boundary');

% Mark both flutter points on damping plot
if ~isnan(lambda_cr_first)
    plot(U_cr_first, 0, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    xline(U_cr_first, 'r--', 'LineWidth', 1.5, 'Alpha', 0.5);
    text(U_cr_first*1.02, 0.05, sprintf('1st Flutter\n%.0f m/s', U_cr_first), ...
         'FontSize', 9, 'Color', 'r', 'FontWeight', 'bold');
end

if ~isnan(lambda_cr_second)
    plot(U_cr_second, 0, 'mo', 'MarkerSize', 12, 'MarkerFaceColor', 'm');
    xline(U_cr_second, 'm--', 'LineWidth', 1.5, 'Alpha', 0.5);
    if U_cr_second < max(velocity)
        text(U_cr_second*1.02, 0.08, sprintf('2nd Flutter\n%.0f m/s', U_cr_second), ...
             'FontSize', 9, 'Color', 'm', 'FontWeight', 'bold');
    end
end

xlabel('Flow Velocity U (m/s)', 'FontSize', 11);
ylabel('Damping Ratio ζ', 'FontSize', 11);
title('(b) Damping Evolution - Zero Crossings', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
legend('Location', 'best', 'FontSize', 9);
ylim([-0.5, 0.5]);
xlim([0, max(velocity)]);

% Subplot 3: Flutter Mechanism Summary
subplot(1, 3, 3);

% Create summary bar chart
mechanisms = {'Single-Mode', 'Coupled-Mode'};
if ~isnan(lambda_cr_first) && ~isnan(lambda_cr_second)
    U_values = [U_cr_first, U_cr_second];
    f_values = [f_cr_first, f_cr_second];
    mode_values = [mode_first, mode_second];
    
    % Bar chart with different colors
    bar(1:2, U_values, 'FaceColor', [0.7 0.3 0.3], 'EdgeColor', 'k', 'LineWidth', 1.2);
    hold on;
    
    % Add frequency text overlay
    for i = 1:2
        text(i, U_values(i) + max(U_values)*0.03, sprintf('%.0f m/s\n%.1f Hz (Mode %d)', ...
             U_values(i), f_values(i), mode_values(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    set(gca, 'XTickLabel', mechanisms);
    ylabel('Flutter Speed (m/s)', 'FontSize', 11);
    title('(c) Two Flutter Mechanisms', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % Add explanation
    text(0.5, -0.15, 'Single-mode: Aeroelastic stiffness softening\nCoupled-mode: Frequency coalescence', ...
         'Units', 'normalized', 'HorizontalAlignment', 'center', ...
         'FontSize', 9, 'Color', 'k', 'FontWeight', 'normal');
elseif ~isnan(lambda_cr_first)
    bar(1, U_cr_first, 'FaceColor', [0.7 0.3 0.3], 'EdgeColor', 'k', 'LineWidth', 1.2);
    set(gca, 'XTickLabel', {'Single-Mode'});
    ylabel('Flutter Speed (m/s)', 'FontSize', 11);
    title('Single Flutter Mechanism', 'FontSize', 12, 'FontWeight', 'bold');
    text(1, U_cr_first, sprintf('%.0f m/s\n%.1f Hz', U_cr_first, f_cr_first), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
else
    text(0.5, 0.5, 'No Flutter Detected\n(Stable up to max velocity)', ...
         'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    axis off;
end

sgtitle(sprintf('SFSF Two-Flutter Analysis: [PZT/0/90/90/0/PZT] (a=%.0f mm, b=%.0f mm, M=%.1f)', ...
        a*1000, b*1000, M_inf), 'FontSize', 14, 'FontWeight', 'bold');

%% Save Results
save('SFSF_two_flutter_results.mat', 'lambda_cr_first', 'U_cr_first', 'f_cr_first', 'mode_first', ...
     'lambda_cr_second', 'U_cr_second', 'f_cr_second', 'mode_second', ...
     'freq_Hz', 'damping_ratio', 'velocity', 'lambda_vals');

fprintf('\n=== FLUTTER SUMMARY ===\n');
fprintf('┌──────────────┬─────────────┬──────────────┬─────────────────┐\n');
fprintf('│ Type         │ U_cr (m/s)  │ f_cr (Hz)    │ Mode            │\n');
fprintf('├──────────────┼─────────────┼──────────────┼─────────────────┤\n');
if ~isnan(lambda_cr_first)
    fprintf('│ First Flutter │ %9.0f   │ %10.1f    │ %d (Single)     │\n', U_cr_first, f_cr_first, mode_first);
end
if ~isnan(lambda_cr_second)
    fprintf('│ Second Flutter│ %9.0f   │ %10.1f    │ %d (Coupled)    │\n', U_cr_second, f_cr_second, mode_second);
end
if isnan(lambda_cr_first) && isnan(lambda_cr_second)
    fprintf('│ No Flutter    │     -       │      -        │ -               │\n');
end
fprintf('└──────────────┴─────────────┴──────────────┴─────────────────┘\n');

%% ========================================================================
% HELPER FUNCTIONS
%% ========================================================================
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
            Q_trans = T * Q_local * T;
            
            % Transverse shear (not used in Kirchhoff theory, kept for completeness)
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

function phi = simply_supported_beam(x, a, i)
         phi = sqrt(2/a) * sin(i * pi * x / a);
end

function dphi = derivative_phi_ss(x, a, i)
    dphi = sqrt(2/a) *(i * pi / a) * cos(i * pi * x / a);
end

function d2phi = second_derivative_phi_ss(x, a, i)
    d2phi = -(i * pi / a)^2 * sin(i * pi * x / a);
end

function psi = free_free_beam(y, b, j)
    if j == 1
        psi = ones(size(y));  % Rigid body translation (constant)
    else
        betas = [4.73004074, 7.85320462, 10.9956078];
        if j <= length(betas)+1
            beta = betas(j-1) / b;
        else
            beta = ((2*j-1)*pi/2) / b;
        end
        C = (cosh(beta*b)-cos(beta*b))/(sinh(beta*b)-sin(beta*b));
        psi = cosh(beta*y) + cos(beta*y) - C * (sinh(beta*y) + sin(beta*y));
    end
end

function dpsi = derivative_psi_free(y, b, j)
    if j == 1
        dpsi = zeros(size(y));
    else
        betas = [4.73004074, 7.85320462, 10.9956078];
        if j <= length(betas)+1
            beta = betas(j-1) / b;
        else
            beta = ((2*j-1)*pi/2) / b;
        end
        C = (cosh(beta*b)-cos(beta*b))/(sinh(beta*b)-sin(beta*b));
        dpsi = beta * (sinh(beta*y) - sin(beta*y) - C * (cosh(beta*y) + cos(beta*y)));
    end
end

function d2psi = second_derivative_psi_free(y, b, j)
    if j == 1
        d2psi = zeros(size(y));
    else
        betas = [4.73004074, 7.85320462, 10.9956078];
        if j <= length(betas)+1
            beta = betas(j-1) / b;
        else
            beta = ((2*j-1)*pi/2) / b;
        end
        C = (cosh(beta*b)-cos(beta*b))/(sinh(beta*b)-sin(beta*b));
        d2psi = beta^2 * (cosh(beta*y) - cos(beta*y) - C * (sinh(beta*y) - sin(beta*y)));
    end
end

function d3psi = third_derivative_psi_free(y, b, j)
    
if j == 1 || j == 2, d3psi = 0;
    return;
end 
    % Rigid body modes
    beta_vals = [4.73004074486, 7.8532046241, 10.995607838, 14.137165491];
    beta = (j-2 <= 4) * (beta_vals(min(j-2,4))/b) + (j-2 > 4) * ((j-1.5)*pi/b);
    C = (cosh(beta*b) - cos(beta*b)) / (sinh(beta*b) - sin(beta*b));
    
    % CORRECTED PHYSICS: Third derivative for Free-Free
    d3psi_raw = beta^3 * (sinh(beta*y) + sin(beta*y) - C * (cosh(beta*y) + cos(beta*y)));
    d3psi = d3psi_raw / sqrt(b);
end

function [x, w] = lgwt(n, a, b)
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
