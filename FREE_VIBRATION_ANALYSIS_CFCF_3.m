%% FREE VIBRATION ANALYSIS OF 6-LAYER COMPOSITE PLATE (CFCF)
% Boundary Conditions: Clamped at x=0, Clamped at x=a, Free at y=0, Free at y=b

clear; clc; close all;

%% 1. PLATE GEOMETRY
a = 0.4;      % Length in x-direction (m) - CLAMPED at x=0 and x=a
b = 0.2;      % Width in y-direction (m) - FREE at y=0 and y=b

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

fprintf('=== LAMINATE PROPERTIES ===\n');
fprintf('D11 = %.2f N·m\n', D11);
fprintf('D22 = %.2f N·m\n', D22);
fprintf('D12 = %.2f N·m\n', D12);
fprintf('D66 = %.2f N·m\n', D66);
fprintf('I0 = %.4f kg/m²\n', I0_mass);
fprintf('Total thickness = %.4f m (%.1f mm)\n', h_total, h_total*1000);
fprintf('\n');

%% 4. RAYLEIGH-RITZ METHOD FOR CFCF
% x-direction: Clamped-Clamped beam functions
% y-direction: Free-Free beam functions (constant for first mode)

% Number of terms
m_terms = 5;  % Terms in x-direction
n_terms = 5;  % Terms in y-direction
total_modes = m_terms * n_terms;

%% ========================================================================
% SECTION: ENHANCED FLUTTER DETECTION FOR CFCF PLATE (USING YOUR FUNCTIONS)
% ========================================================================
% This code uses your corrected beam functions:
%   - clamped_clamped_beam(x, a, i)
%   - derivative_phi_corrected(x, a, i)
%   - second_derivative_phi_corrected(x, a, i)
%   - free_free_beam_corrected(y, b, j)
%   - derivative_psi_corrected(y, b, j)
%   - second_derivative_psi_corrected(y, b, j)

fprintf('\n========================================\n');
fprintf('ENHANCED FLUTTER DETECTION: CFCF PLATE\n');
fprintf('Using Corrected Beam Functions\n');
fprintf('========================================\n');

% Check if required variables exist
if ~exist('M_modal', 'var') || ~exist('K_modal', 'var')
    error('Required variables not defined. Run main flutter analysis first.');
end

% Extended λ range to capture second flutter (CFCF needs higher λ)
lambda_min = 0;
lambda_max = 3500;  % CFCF is stiffer, needs higher λ for 2nd flutter
n_lambda = 800;     % Number of points
lambda_values_enh = linspace(lambda_min, lambda_max, n_lambda);

% Storage for all 4 modes
freq_Hz_all = zeros(4, n_lambda);
damping_all = zeros(4, n_lambda);

% Flutter detection structures
flutter_pairs = {};
flutter_count = 0;

fprintf('CFCF Configuration:\n');
fprintf('  - Clamped at x=0 and x=a (both ends)\n');
fprintf('  - Free at y=0 and y=b\n');
fprintf('  - λ range: 0 to %.0f (%d points)\n', lambda_max, n_lambda);
fprintf('  - Expected 1st flutter: λ ≈ 330-350\n');
fprintf('  - Expected 2nd flutter: λ ≈ 2500-2700\n');

fprintf('\nComputing eigenvalue evolution for CFCF plate...\n');

for k = 1:n_lambda
    lambda_val = lambda_values_enh(k);
    
    % Dynamic pressure
    q_dyn_val = lambda_val * beta_flow * D11 / (2 * a^3);
    U_val = sqrt(2 * q_dyn_val / rho_air);
    
    % Aerodynamic matrices (4x4)
    K_aero_val = (2 * q_dyn_val / beta_flow) * A_aero;
    
    % Aerodynamic damping (piston theory)
    if U_val > 0 && M_inf^2 > 2
        g_a_val = (rho_air * U_val * (M_inf^2 - 2)) / (I0_mass * beta_flow^3);
        C_aero_val = g_a_val * M_modal;
    else
        C_aero_val = zeros(4);
    end
    
    % Small structural damping (0.5%)
    zeta_struct = 0.005;
    C_struct_val = 2 * zeta_struct * diag(omega_n) * M_modal;
    
    % Total matrices
    M_sys = M_modal;
    K_sys = K_modal + K_aero_val;
    C_sys = C_struct_val + C_aero_val;
    
    % State-space eigenvalue problem (8x8)
    A_state = [zeros(4), eye(4); -M_sys\K_sys, -M_sys\C_sys];
    [~, D_state] = eig(A_state);
    eigenvalues = diag(D_state);
    
    % Extract oscillatory modes
    eig_osc = eigenvalues(imag(eigenvalues) > 1e-6 & abs(real(eigenvalues)) < abs(imag(eigenvalues)));
    
    if ~isempty(eig_osc)
        [~, sort_idx] = sort(imag(eig_osc));
        eig_osc = eig_osc(sort_idx);
        
        for j = 1:min(4, length(eig_osc))
            freq_Hz_all(j, k) = imag(eig_osc(j)) / (2*pi);
            damping_all(j, k) = -real(eig_osc(j)) / abs(eig_osc(j));
        end
    end
    
    % --- Flutter Detection for Mode Pairs ---
    % Pair 1: Modes 1-2 (classical bending-torsion)
    if k > 10
        freq_diff_12 = abs(freq_Hz_all(2, k) - freq_Hz_all(1, k));
        if freq_diff_12 < 3.0 && freq_diff_12 > 0 && ~any(strcmp(flutter_pairs, 'Pair1-2'))
            flutter_count = flutter_count + 1;
            flutter_pairs{flutter_count} = 'Pair1-2';
            lambda_cr_12 = lambda_values_enh(k);
            freq_flutter_12 = (freq_Hz_all(1, k) + freq_Hz_all(2, k)) / 2;
            
            % Calculate velocity and Mach
            q_dyn_cr = lambda_cr_12 * beta_flow * D11 / (2 * a^3);
            U_cr_12 = sqrt(2 * q_dyn_cr / rho_air);
            Mach_cr_12 = U_cr_12 / 340;
            
            fprintf('\n✓ FIRST FLUTTER DETECTED (Modes 1-2):\n');
            fprintf('  λ_cr = %.1f\n', lambda_cr_12);
            fprintf('  Flutter velocity = %.0f m/s (Mach %.2f)\n', U_cr_12, Mach_cr_12);
            fprintf('  Flutter frequency = %.1f Hz\n', freq_flutter_12);
        end
    end
    
    % Pair 2: Modes 3-4 (higher mode flutter)
    if k > 10 && freq_Hz_all(3, k) > 0 && freq_Hz_all(4, k) > 0
        freq_diff_34 = abs(freq_Hz_all(4, k) - freq_Hz_all(3, k));
        if freq_diff_34 < 3.0 && freq_diff_34 > 0 && ~any(strcmp(flutter_pairs, 'Pair3-4'))
            flutter_count = flutter_count + 1;
            flutter_pairs{flutter_count} = 'Pair3-4';
            lambda_cr_34 = lambda_values_enh(k);
            freq_flutter_34 = (freq_Hz_all(3, k) + freq_Hz_all(4, k)) / 2;
            
            % Calculate velocity and Mach
            q_dyn_cr_34 = lambda_cr_34 * beta_flow * D11 / (2 * a^3);
            U_cr_34 = sqrt(2 * q_dyn_cr_34 / rho_air);
            Mach_cr_34 = U_cr_34 / 340;
            
            fprintf('\n✓ SECOND FLUTTER DETECTED (Modes 3-4):\n');
            fprintf('  λ_cr = %.1f\n', lambda_cr_34);
            fprintf('  Flutter velocity = %.0f m/s (Mach %.2f)\n', U_cr_34, Mach_cr_34);
            fprintf('  Flutter frequency = %.1f Hz\n', freq_flutter_34);
        end
    end
    
    if mod(k, 100) == 0
        fprintf('  Progress: %.0f%% (λ = %.0f)\n', k/n_lambda*100, lambda_val);
    end
end

%% ========================================================================
% SECTION: ENHANCED PLOT FOR CFCF (Both Flutter Points)
% ========================================================================
figure('Color', 'w', 'Position', [100, 100, 1600, 700]);

% Subplot 1: Frequency Coalescence (All 4 modes)
subplot(1,2,1);
colors = {'b', 'r', 'g', 'm'};
mode_names = {'Mode 1 (1st Bending, 22.9 Hz)', ...
              'Mode 2 (1st Torsion, 167.6 Hz)', ...
              'Mode 3 (2nd Bending, 539.4 Hz)', ...
              'Mode 4 (2nd Torsion, 615.6 Hz)'};

for j = 1:4
    valid_idx = freq_Hz_all(j, :) > 0;
    plot(lambda_values_enh(valid_idx), freq_Hz_all(j, valid_idx), colors{j}, ...
         'LineWidth', 2, 'DisplayName', mode_names{j});
    hold on;
end

% Mark first flutter (Modes 1-2)
if exist('lambda_cr_12', 'var')
    plot(lambda_cr_12, freq_flutter_12, 'ko', 'MarkerSize', 14, ...
         'MarkerFaceColor', 'k', 'LineWidth', 2, ...
         'DisplayName', sprintf('1st Flutter: λ=%.0f, f=%.0fHz', lambda_cr_12, freq_flutter_12));
    xline(lambda_cr_12, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

% Mark second flutter (Modes 3-4)
if exist('lambda_cr_34', 'var')
    plot(lambda_cr_34, freq_flutter_34, 'ro', 'MarkerSize', 14, ...
         'MarkerFaceColor', 'r', 'LineWidth', 2, ...
         'DisplayName', sprintf('2nd Flutter: λ=%.0f, f=%.0fHz', lambda_cr_34, freq_flutter_34));
    xline(lambda_cr_34, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

xlabel('\lambda (Aerodynamic Pressure Parameter)', 'FontSize', 12);
ylabel('Frequency (Hz)', 'FontSize', 12);
title('(a) CFCF: Frequency Coalescence - Multiple Flutter Mechanisms', ...
    'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, min(lambda_max, 3000)]);
ylim([0, 800]);

% Subplot 2: Damping Evolution
subplot(1,2,2);
for j = 1:4
    valid_idx = damping_all(j, :) ~= 0;
    plot(lambda_values_enh(valid_idx), damping_all(j, valid_idx)*100, colors{j}, ...
         'LineWidth', 2, 'DisplayName', mode_names{j});
    hold on;
end

% Mark flutter points on damping plot
if exist('lambda_cr_12', 'var')
    plot(lambda_cr_12, 0, 'ko', 'MarkerSize', 14, 'MarkerFaceColor', 'k', 'LineWidth', 2);
    xline(lambda_cr_12, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end
if exist('lambda_cr_34', 'var')
    plot(lambda_cr_34, 0, 'ro', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'LineWidth', 2);
    xline(lambda_cr_34, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('\lambda', 'FontSize', 12);
ylabel('Damping Ratio ζ (%)', 'FontSize', 12);
title('(b) CFCF: Damping Evolution - Multiple Zero Crossings', ...
    'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, min(lambda_max, 3000)]);
ylim([-2, 15]);

sgtitle('CFCF Plate: Multiple Flutter Mechanisms (Modes 1-2 and 3-4)', ...
    'FontSize', 14, 'FontWeight', 'bold');

%% ========================================================================
% SECTION: SAVE FIGURE
% ========================================================================
fprintf('\n========================================\n');
fprintf('SAVING FIGURE\n');
fprintf('========================================\n');

try
    exportgraphics(gcf, 'CFCF_Enhanced_Flutter_Analysis.png', 'Resolution', 300);
    fprintf('✓ Figure saved: CFCF_Enhanced_Flutter_Analysis.png\n');
catch
    try
        print(gcf, 'CFCF_Enhanced_Flutter_Analysis.png', '-dpng', '-r300');
        fprintf('✓ Figure saved: CFCF_Enhanced_Flutter_Analysis.png (using print)\n');
    catch
        fprintf('⚠ Auto-save failed. Please save manually:\n');
        fprintf('   File → Save As → PNG → "CFCF_Enhanced_Flutter_Analysis.png"\n');
    end
end

%% ========================================================================
% SECTION: SUMMARY TABLE
% ========================================================================
fprintf('\n========================================\n');
fprintf('CFCF FLUTTER ANALYSIS - COMPLETE SUMMARY\n');
fprintf('========================================\n');

fprintf('\n┌─────────────────────────────────────────────────────────────────┐\n');
fprintf('│                    CFCF PLATE PROPERTIES                        │\n');
fprintf('├─────────────────────────────────────────────────────────────────┤\n');
fprintf('│ Length a (x-direction):        %.3f m (%.0f mm)                 │\n', a, a*1000);
fprintf('│ Width b (y-direction):         %.3f m (%.0f mm)                 │\n', b, b*1000);
fprintf('│ Aspect ratio a/b:              %.2f                            │\n', a/b);
fprintf('│ D11 (bending stiffness):       %.2f N·m                        │\n', D11);
fprintf('│ I0 (mass per area):            %.4f kg/m²                      │\n', I0_mass);
fprintf('└─────────────────────────────────────────────────────────────────┘\n');

fprintf('\n┌─────────────────────────────────────────────────────────────────┐\n');
fprintf('│                    NATURAL FREQUENCIES                          │\n');
fprintf('├─────────────────────────────────────────────────────────────────┤\n');
for i = 1:4
    fprintf('│ Mode %d:                         %.1f Hz                         │\n', i, freq_n(i));
end
fprintf('└─────────────────────────────────────────────────────────────────┘\n');

fprintf('\n┌─────────────────────────────────────────────────────────────────┐\n');
fprintf('│                    FLUTTER RESULTS                              │\n');
fprintf('├─────────────────────────────────────────────────────────────────┤\n');
if exist('lambda_cr_12', 'var')
    fprintf('│ 1st Flutter (Modes 1-2):                                      │\n');
    fprintf('│   λ_cr = %.1f                                                │\n', lambda_cr_12);
    fprintf('│   U_cr = %.0f m/s (Mach %.2f)                                │\n', U_cr_12, Mach_cr_12);
    fprintf('│   f_cr = %.1f Hz                                             │\n', freq_flutter_12);
end
if exist('lambda_cr_34', 'var')
    fprintf('│ 2nd Flutter (Modes 3-4):                                      │\n');
    fprintf('│   λ_cr = %.1f                                                │\n', lambda_cr_34);
    fprintf('│   U_cr = %.0f m/s (Mach %.2f)                                │\n', U_cr_34, Mach_cr_34);
    fprintf('│   f_cr = %.1f Hz                                             │\n', freq_flutter_34);
else
    fprintf('│ 2nd Flutter (Modes 3-4):   Not detected in λ range           │\n');
end
fprintf('└─────────────────────────────────────────────────────────────────┘\n');

fprintf('\n========================================\n');
fprintf('ENHANCED FLUTTER ANALYSIS COMPLETE\n');
fprintf('========================================\n');

%% Helper functions (include all necessary functions at the end)
% (Include calculate_laminate_properties, lgwt, beam functions, etc.)
%% ========================================================================
% SECTION 12: LAMINATE PROPERTIES CALCULATION
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
            Q_trans = T * Q_local * T;
            
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
function [x, w] = lgwt(N, a, b)
    % Legendre-Gauss quadrature nodes and weights
    N = N-1;
    N1 = N+1;
    N2 = N+2;
    
    xu = linspace(-1, 1, N1)';
    
    % Initial guess
    y = cos((2*(0:N)'+1)*pi/(2*N+2)) + (0.27/N1)*sin(pi*xu*N/N2);
    
    % Legendre-Gauss Vandermonde matrix
    L = zeros(N1, N2);
    
    % Derivative of LG
    Lp = zeros(N1, N2);
    
    % Loop for iteration
    for i = 1:N1
        for j = 1:N1
            if j == 1
                L(i,j) = 1;
                Lp(i,j) = 0;
            elseif j == 2
                L(i,j) = y(i);
                Lp(i,j) = 1;
            else
                L(i,j) = ((2*j-3)*y(i)*L(i,j-1) - (j-2)*L(i,j-2))/(j-1);
                Lp(i,j) = ((2*j-3)*(L(i,j-1) + y(i)*Lp(i,j-1)) - (j-2)*Lp(i,j-2))/(j-1);
            end
        end
    end
    
    x = (a*(1-y) + b*(1+y))/2;
    w = (b-a) ./ ((1-y.^2) .* Lp(:,N1).^2) * (N2/N1)^2;
end

%% ========================================================================
% CORRECT BEAM FUNCTIONS - SINGLE CONSISTENT VERSION
% ========================================================================

function phi = clamped_clamped_beam(x, a, i)
    % Clamped-Clamped beam - AUTO-NORMALIZED for any a
    beta_vals = [4.73004074486, 7.8532046241, 10.995607838, 14.137165491, ...
                 17.278759657, 20.420352245, 23.561944902, 26.703537555];
    
    if i <= length(beta_vals)
        beta = beta_vals(i) / a;
    else
        beta = (i + 0.5) * pi / a;
    end
    
    C = (cosh(beta*a) - cos(beta*a)) / (sinh(beta*a) - sin(beta*a));
    phi_raw = cosh(beta*x) - cos(beta*x) - C * (sinh(beta*x) - sin(beta*x));
    
    % Auto-compute normalization
    persistent norm_cc
    if isempty(norm_cc)
        norm_cc = zeros(8,1);
        for ii = 1:8
            if ii <= length(beta_vals)
                b_i = beta_vals(ii) / a;
            else
                b_i = (ii + 0.5) * pi / a;
            end
            C_i = (cosh(b_i*a) - cos(b_i*a)) / (sinh(b_i*a) - sin(b_i*a));
            f = @(xx) (cosh(b_i*xx) - cos(b_i*xx) - C_i*(sinh(b_i*xx) - sin(b_i*xx))).^2;
            norm_cc(ii) = sqrt(integral(f, 0, a));
        end
    end
    
    if i <= length(norm_cc)
        phi = phi_raw / norm_cc(i);
    else
        phi = phi_raw;
    end
end

function dphi = derivative_phi_corrected(x, a, i)
    % First derivative - uses same normalization as base function
    beta_vals = [4.73004074486, 7.8532046241, 10.995607838, 14.137165491];
    
    if i <= length(beta_vals)
        beta = beta_vals(i) / a;
    else
        beta = (i + 0.5) * pi / a;
    end
    
    C = (cosh(beta*a) - cos(beta*a)) / (sinh(beta*a) - sin(beta*a));
    dphi_raw = beta * sinh(beta*x) + beta * sin(beta*x) - ...
               C * beta * (cosh(beta*x) - cos(beta*x));
    
    % Get same normalization as base function
    persistent norm_cc
    if isempty(norm_cc)
        norm_cc = zeros(8,1);
        for ii = 1:8
            if ii <= length(beta_vals)
                b_i = beta_vals(ii) / a;
            else
                b_i = (ii + 0.5) * pi / a;
            end
            C_i = (cosh(b_i*a) - cos(b_i*a)) / (sinh(b_i*a) - sin(b_i*a));
            f = @(xx) (cosh(b_i*xx) - cos(b_i*xx) - C_i*(sinh(b_i*xx) - sin(b_i*xx))).^2;
            norm_cc(ii) = sqrt(integral(f, 0, a));
        end
    end
    
    if i <= length(norm_cc)
        dphi = dphi_raw / norm_cc(i);
    else
        dphi = dphi_raw;
    end
end

function d2phi = second_derivative_phi_corrected(x, a, i)
    % Second derivative - uses same normalization as base function
    beta_vals = [4.73004074486, 7.8532046241, 10.995607838, 14.137165491];
    
    if i <= length(beta_vals)
        beta = beta_vals(i) / a;
    else
        beta = (i + 0.5) * pi / a;
    end
    
    C = (cosh(beta*a) - cos(beta*a)) / (sinh(beta*a) - sin(beta*a));
    d2phi_raw = beta^2 * cosh(beta*x) + beta^2 * cos(beta*x) - ...
                C * beta^2 * (sinh(beta*x) + sin(beta*x));
    
    persistent norm_cc
    if isempty(norm_cc)
        norm_cc = zeros(8,1);
        for ii = 1:8
            if ii <= length(beta_vals)
                b_i = beta_vals(ii) / a;
            else
                b_i = (ii + 0.5) * pi / a;
            end
            C_i = (cosh(b_i*a) - cos(b_i*a)) / (sinh(b_i*a) - sin(b_i*a));
            f = @(xx) (cosh(b_i*xx) - cos(b_i*xx) - C_i*(sinh(b_i*xx) - sin(b_i*xx))).^2;
            norm_cc(ii) = sqrt(integral(f, 0, a));
        end
    end
    
    if i <= length(norm_cc)
        d2phi = d2phi_raw / norm_cc(i);
    else
        d2phi = d2phi_raw;
    end
end

%% CORRECTED FREE-FREE BEAM FUNCTIONS
function psi = free_free_beam_corrected(y, b, j)
    if j == 1
        psi = 1 / sqrt(b); % Rigid body translation, unit norm
        return;
    end
    
    beta_vals = [4.73004074486, 7.8532046241, 10.995607838, 14.137165491];
    if j-1 <= length(beta_vals)
        beta = beta_vals(j-1) / b;
    else
        beta = (j - 0.5) * pi / b;
    end
    C = (cosh(beta*b) - cos(beta*b)) / (sinh(beta*b) - sin(beta*b));
    psi_raw = cosh(beta*y) + cos(beta*y) - C * (sinh(beta*y) + sin(beta*y));
    psi = psi_raw / sqrt(b); % ∫ψ²dy = b → normalization by √b gives unit mass
end

function dpsi = derivative_psi_corrected(y, b, j)
    if j == 1, dpsi = 0; return; end
    beta_vals = [4.73004074486, 7.8532046241, 10.995607838, 14.137165491];
    if j-1 <= length(beta_vals)
        beta = beta_vals(j-1) / b;
    else
        beta = (j - 0.5) * pi / b;
    end
    C = (cosh(beta*b) - cos(beta*b)) / (sinh(beta*b) - sin(beta*b));
    
    % CORRECTED: -C*(cosh + cos) instead of -C*(cosh - cos)
    dpsi_raw = beta * (sinh(beta*y) - sin(beta*y) - C * (cosh(beta*y) + cos(beta*y)));
    dpsi = dpsi_raw / sqrt(b); % Same normalization as base function
end

function d2psi = second_derivative_psi_corrected(y, b, j)
    if j == 1, d2psi = 0; return; end
    beta_vals = [4.73004074486, 7.8532046241, 10.995607838, 14.137165491];
    if j-1 <= length(beta_vals)
        beta = beta_vals(j-1) / b;
    else
        beta = (j - 0.5) * pi / b;
    end
    C = (cosh(beta*b) - cos(beta*b)) / (sinh(beta*b) - sin(beta*b));
    
    d2psi_raw = beta^2 * (cosh(beta*y) - cos(beta*y) - C * (sinh(beta*y) - sin(beta*y)));
    d2psi = d2psi_raw / sqrt(b);
end

function d3psi = third_derivative_psi_corrected(y, b, j)
    if j == 1, d3psi = 0; return; end
    beta_vals = [4.73004074486, 7.8532046241, 10.995607838, 14.137165491];
    if j-1 <= length(beta_vals)
        beta = beta_vals(j-1) / b;
    else
        beta = (j - 0.5) * pi / b;
    end
    C = (cosh(beta*b) - cos(beta*b)) / (sinh(beta*b) - sin(beta*b));
    
    % CORRECTED: -C*(cosh - cos) instead of -C*(cosh + cos)
    d3psi_raw = beta^3 * (sinh(beta*y) + sin(beta*y) - C * (cosh(beta*y) - cos(beta*y)));
    d3psi = d3psi_raw / sqrt(b);
end
