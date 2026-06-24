%% Flutter Analysis for PZT/0/90/90/0/PZT Laminate
% Clear workspace
clear; clc; close all;

%% Define geometry and material properties
a = 0.3;           % Length (m)
b = 0.2;           % Width (m)
alpha = a/b;          % Aspect ratio (approximately 0.292)

% Reference values for normalization
E_ref = 1e9;          % 1 GPa reference
rho_ref = 1000;       % 1000 kg/m³ reference
eps_ref = 8.854e-12;  % Permittivity of free space (F/m)

%% Material definitions
% PZT-5A layer (isotropic, piezoelectric)
pzt = struct('type', 'piezo', ...
             'E', 66e9, ...           % Young's modulus (Pa)
             'nu', 0.31, ...          % Poisson's ratio
             'rho', 7800, ...         % Density (kg/m³)
             'thickness', 0.0005, ... % 0.5 mm thickness
             'e31', -5.4, ...         % Piezoelectric constant (C/m²)
             'eps33', 1.5e-8);        % Permittivity (F/m)

% Composite layer (0° orientation)
comp0 = struct('type', 'comp', ...
               'E1', 140e9, ...       % Longitudinal modulus (Pa)
               'E2', 10e9, ...        % Transverse modulus (Pa)
               'G12', 5e9, ...        % Shear modulus (Pa)
               'nu12', 0.3, ...       % Poisson's ratio
               'rho', 1600, ...       % Density (kg/m³)
               'thickness', 0.001, ...% 1 mm thickness
               'angle', 0);           % 0° orientation

% Composite layer (90° orientation)
comp90 = comp0;
comp90.angle = 90;

% Laminate stacking: [PZT/0/90/90/0/PZT] (symmetric)
layers = {pzt, comp0, comp90, comp90, comp0, pzt};

%% Calculate laminate properties
[A, B, D, As, I0, I2] = calculate_laminate_properties(layers, a, b);

% Extract relevant stiffness components
D11 = D(1,1);          % Bending stiffness in x-direction
D22 = D(2,2);          % Bending stiffness in y-direction
D12 = D(1,2);          % Poisson coupling
D66 = D(3,3);          % Twisting stiffness

% Shear stiffnesses (from As matrix)
A55_s = As(1,1);       % Transverse shear stiffness (x-z plane)
A44_s = As(2,2);       % Transverse shear stiffness (y-z plane)

% Mass moments of inertia
I0_mass = I0;          % Mass per unit area (kg/m²)
I2_mass = I2;          % Mass moment of inertia (kg)

%% Normalize parameters (as per Eq. ~)
E_ref_norm = 1e9;      % 1 GPa
rho_ref_norm = 1000;   % 1000 kg/m³

% Normalized material properties (for reference)
E1_norm = comp0.E1 / E_ref_norm;
E2_norm = comp0.E2 / E_ref_norm;
G12_norm = comp0.G12 / E_ref_norm;
rho_norm = comp0.rho / rho_ref_norm;
e31_norm = pzt.e31 / sqrt(E_ref_norm * eps_ref);
eps33_norm = pzt.eps33 / eps_ref;

%% Dimensionless system parameters (as per Section 3.2)
% Reference quantities for normalization
h_total = sum(cellfun(@(x) x.thickness, layers));  % Total thickness
omega_ref = sqrt(D11 / (I0_mass * a^4));           % Reference frequency

% Normalized parameters
D11_bar = 1.0;  % Reference bending stiffness
D22_bar = D22 / D11;
D12_bar = D12 / D11;
D66_bar = D66 / D11;

A55_s_bar = A55_s * a^2 / D11;
A44_s_bar = A44_s * a^2 / D11;

I0_bar = I0_mass * a^2 * omega_ref^2 / D11;
I2_bar = I2_mass * omega_ref^2 / D11;

%% Initial in-plane loads (from piezoelectric actuation)
% For V = 0 (baseline case)
V_applied = 0;  % Applied voltage (Volts)

% Piezoelectric induced forces
% For symmetric PZT layers with same polarization
Nx_pzt = 0; Ny_pzt = 0;  % Initialize
for i = 1:length(layers)
    if strcmp(layers{i}.type, 'piezo')
        % For a plate with clamped edges, induced forces depend on voltage
        % Simplified: N_x^p = N_y^p = e31 * V (for biaxial stress state)
        Nx_pzt = Nx_pzt + layers{i}.e31 * V_applied;
        Ny_pzt = Ny_pzt + layers{i}.e31 * V_applied;
    end
end

% Dimensionless in-plane loads
Nx0_bar = Nx_pzt * a^2 / D11;
Ny0_bar = Ny_pzt * a^2 / D11;
Nxy0_bar = 0;  % No shear load

%% Flow parameters
M = 2.0;                    % Mach number
beta = sqrt(M^2 - 1);       % Prandtl-Glauert factor
U_inf = 111490;             % Flutter velocity (m/s) from analysis
rho_air = 1.225;            % Air density (kg/m³) at sea level

% Dynamic pressure parameter lambda
% From the analysis: lambda_cr ≈ 494.3 at flutter
lambda_cr = 494.3;          % Critical aerodynamic pressure parameter

% Flow speed parameter (dimensionless)
U_bar = U_inf / (a * omega_ref);

% Aerodynamic pressure parameter
lambda_values = linspace(0, 550, 150);  % Range of lambda values

%% Define stiffness matrix coefficients (Mode 1)
pi_val = pi;

% Mode 1 coefficients (from Appendix B, with proper normalization)
k11_1 = A55_s_bar * (pi_val/alpha)^2 + A44_s_bar * pi_val^2 + ...
        Nx0_bar * (pi_val/alpha)^2 + 2*Nxy0_bar * (pi_val^2/alpha) + ...
        Ny0_bar * pi_val^2 - lambda_cr * (pi_val/alpha);

k12_1 = A55_s_bar * (pi_val/alpha);
k13_1 = A44_s_bar * pi_val;
k21_1 = -A55_s_bar * (pi_val/alpha);
k22_1 = D11_bar * (pi_val/alpha)^2 + D66_bar * pi_val^2 + A55_s_bar;
k23_1 = (D12_bar + D66_bar) * (pi_val^2/alpha);
k31_1 = -A44_s_bar * pi_val;
k32_1 = (D12_bar + D66_bar) * (pi_val^2/alpha);
k33_1 = D66_bar * (pi_val/alpha)^2 + D22_bar * pi_val^2 + A44_s_bar;

% Mode 2 coefficients (with factor 2 for second mode shape)
k44_2 = A55_s_bar * (pi_val/alpha)^2 + 4*A44_s_bar * pi_val^2 + ...
        Nx0_bar * (pi_val/alpha)^2 + 4*Nxy0_bar * (pi_val^2/alpha) + ...
        4*Ny0_bar * pi_val^2 - lambda_cr * (pi_val/alpha);

k45_2 = A55_s_bar * (pi_val/alpha);
k46_2 = 2 * A44_s_bar * pi_val;
k54_2 = -A55_s_bar * (pi_val/alpha);
k55_2 = D11_bar * (pi_val/alpha)^2 + 4*D66_bar * pi_val^2 + A55_s_bar;
k56_2 = 2 * (D12_bar + D66_bar) * (pi_val^2/alpha);
k64_2 = -2 * A44_s_bar * pi_val;
k65_2 = 2 * (D12_bar + D66_bar) * (pi_val^2/alpha);
k66_2 = D66_bar * (pi_val/alpha)^2 + 4*D22_bar * pi_val^2 + A44_s_bar;

%% Mass matrix (dimensionless)
M_bar = diag([I0_bar, I2_bar, I2_bar, I0_bar, I2_bar, I2_bar]);

%% CASE 1: With aerodynamic damping
fprintf('Case 1: With aerodynamic damping\n');
fprintf('=================================\n');
fprintf('Laminate properties:\n');
fprintf('  Total thickness: %.2f mm\n', h_total*1000);
fprintf('  D11 = %.2f N·m\n', D11);
fprintf('  I0 = %.4f kg/m²\n', I0_mass);
fprintf('  Reference frequency: %.1f Hz\n', omega_ref/(2*pi));
fprintf('  Aspect ratio (a/b): %.3f\n', alpha);
fprintf('\n');

% Pre-allocate storage
n_lambda = length(lambda_values);
frequencies_damp = zeros(2, n_lambda);
damping_ratios_damp = zeros(2, n_lambda);
freq_Hz_damp = zeros(2, n_lambda);

for i = 1:n_lambda
    lambda = lambda_values(i);
    
    % Aerodynamic pressure effect on stiffness
    lambda_term = lambda * (pi_val/alpha);
    
    % Stiffness matrix (updated for current lambda)
    K_bar = zeros(6,6);
    
    % Mode 1 block
    K_bar(1,1) = k11_1 - lambda_term + Nx0_bar*(pi_val/alpha)^2;
    K_bar(1,2) = k12_1; 
    K_bar(1,3) = k13_1;
    K_bar(2,1) = k21_1; 
    K_bar(2,2) = k22_1; 
    K_bar(2,3) = k23_1;
    K_bar(3,1) = k31_1; 
    K_bar(3,2) = k32_1; 
    K_bar(3,3) = k33_1;
    
    % Mode 2 block
    K_bar(4,4) = k44_2 - lambda_term + 4*Ny0_bar*pi_val^2;
    K_bar(4,5) = k45_2; 
    K_bar(4,6) = k46_2;
    K_bar(5,4) = k54_2; 
    K_bar(5,5) = k55_2; 
    K_bar(5,6) = k56_2;
    K_bar(6,4) = k64_2; 
    K_bar(6,5) = k65_2; 
    K_bar(6,6) = k66_2;
    
    % Damping matrix (aerodynamic damping)
    U_bar_val = 1.0;  % Normalized flow velocity
    C_bar = lambda * U_bar_val * diag([1, 0, 0, 1, 0, 0]);
    
    % State-space form: A * x = iΩ * x
    n = size(M_bar, 1);
    A = [zeros(n), eye(n);
         -M_bar\K_bar, -M_bar\C_bar];
    
    % Solve eigenvalue problem
    [V, D] = eig(A);
    eigenvalues = diag(D);
    
    % Extract frequencies and damping
    % For oscillatory modes, eigenvalues are complex conjugates
    oscillatory_modes = imag(eigenvalues) > 1e-6 & abs(real(eigenvalues)) < 100;
    eig_osc = eigenvalues(oscillatory_modes);
    
    % Sort by frequency
    [~, idx] = sort(abs(imag(eig_osc)));
    
    % Store two lowest frequency modes
    for j = 1:min(2, length(idx))
        mode_idx = idx(j);
        omega_complex = eig_osc(mode_idx);
        
        % Physical frequency (rad/s)
        omega_phys = imag(omega_complex) * omega_ref;
        freq_Hz_damp(j, i) = omega_phys / (2*pi);
        
        % Damping ratio (positive = stable)
        sigma = real(omega_complex);
        damping_ratios_damp(j, i) = -sigma / sqrt(sigma^2 + imag(omega_complex)^2);
    end
end

%% CASE 2: Without aerodynamic damping
fprintf('Case 2: Without aerodynamic damping\n');
fprintf('===================================\n');

frequencies_no_damp = zeros(2, n_lambda);
damping_ratios_no_damp = zeros(2, n_lambda);
freq_Hz_no_damp = zeros(2, n_lambda);

for i = 1:n_lambda
    lambda = lambda_values(i);
    lambda_term = lambda * (pi_val/alpha);
    
    % Stiffness matrix (same as before)
    K_bar = zeros(6,6);
    K_bar(1,1) = k11_1 - lambda_term + Nx0_bar*(pi_val/alpha)^2;
    K_bar(1,2) = k12_1; K_bar(1,3) = k13_1;
    K_bar(2,1) = k21_1; K_bar(2,2) = k22_1; K_bar(2,3) = k23_1;
    K_bar(3,1) = k31_1; K_bar(3,2) = k32_1; K_bar(3,3) = k33_1;
    K_bar(4,4) = k44_2 - lambda_term + 4*Ny0_bar*pi_val^2;
    K_bar(4,5) = k45_2; K_bar(4,6) = k46_2;
    K_bar(5,4) = k54_2; K_bar(5,5) = k55_2; K_bar(5,6) = k56_2;
    K_bar(6,4) = k64_2; K_bar(6,5) = k65_2; K_bar(6,6) = k66_2;
    
    % No damping
    C_bar = zeros(6,6);
    
    % State-space form
    n = size(M_bar, 1);
    A = [zeros(n), eye(n);
         -M_bar\K_bar, -M_bar\C_bar];
    
    [V, D] = eig(A);
    eigenvalues = diag(D);
    
    oscillatory_modes = imag(eigenvalues) > 1e-6 & abs(real(eigenvalues)) < 100;
    eig_osc = eigenvalues(oscillatory_modes);
    
    [~, idx] = sort(abs(imag(eig_osc)));
    
    for j = 1:min(2, length(idx))
        mode_idx = idx(j);
        omega_complex = eig_osc(mode_idx);
        omega_phys = imag(omega_complex) * omega_ref;
        freq_Hz_no_damp(j, i) = omega_phys / (2*pi);
        
        sigma = real(omega_complex);
        damping_ratios_no_damp(j, i) = -sigma / sqrt(sigma^2 + imag(omega_complex)^2);
    end
end

%% Plot results
figure('Position', [50, 50, 1400, 900]);

% Plot 1: Frequency evolution (with damping)
subplot(2,2,1);
plot(lambda_values, freq_Hz_damp(1,:)/1000, 'b-', 'LineWidth', 2); hold on;
plot(lambda_values, freq_Hz_damp(2,:)/1000, 'r-', 'LineWidth', 2);
xlabel('\lambda (Aerodynamic Pressure Parameter)', 'FontSize', 12);
ylabel('Frequency (kHz)', 'FontSize', 12);
title('Frequency Evolution (With Damping)', 'FontSize', 14);
legend('Mode 1', 'Mode 2', 'Location', 'best');
grid on;
xline(lambda_cr, 'k--', 'LineWidth', 1.5);
text(lambda_cr, max(freq_Hz_damp(:))/2000, sprintf('\\lambda_{cr}=%.1f', lambda_cr), ...
     'HorizontalAlignment', 'center', 'FontSize', 10);

% Plot 2: Damping evolution (with damping)
subplot(2,2,2);
plot(lambda_values, damping_ratios_damp(1,:), 'b-', 'LineWidth', 2); hold on;
plot(lambda_values, damping_ratios_damp(2,:), 'r-', 'LineWidth', 2);
xlabel('\lambda', 'FontSize', 12);
ylabel('Damping Ratio', 'FontSize', 12);
title('Damping Evolution (With Damping)', 'FontSize', 14);
legend('Mode 1', 'Mode 2', 'Location', 'best');
grid on;
yline(0, 'k--', 'LineWidth', 1);
xline(lambda_cr, 'k--', 'LineWidth', 1.5);

% Plot 3: Frequency evolution (without damping)
subplot(2,2,3);
plot(lambda_values, freq_Hz_no_damp(1,:)/1000, 'b-', 'LineWidth', 2); hold on;
plot(lambda_values, freq_Hz_no_damp(2,:)/1000, 'r-', 'LineWidth', 2);
xlabel('\lambda', 'FontSize', 12);
ylabel('Frequency (kHz)', 'FontSize', 12);
title('Frequency Evolution (Without Damping)', 'FontSize', 14);
legend('Mode 1', 'Mode 2', 'Location', 'best');
grid on;
xline(lambda_cr, 'k--', 'LineWidth', 1.5);

% Plot 4: Damping evolution (without damping)
subplot(2,2,4);
plot(lambda_values, damping_ratios_no_damp(1,:), 'b-', 'LineWidth', 2); hold on;
plot(lambda_values, damping_ratios_no_damp(2,:), 'r-', 'LineWidth', 2);
xlabel('\lambda', 'FontSize', 12);
ylabel('Damping Ratio', 'FontSize', 12);
title('Damping Evolution (Without Damping)', 'FontSize', 14);
legend('Mode 1', 'Mode 2', 'Location', 'best');
grid on;
yline(0, 'k--', 'LineWidth', 1);
xline(lambda_cr, 'k--', 'LineWidth', 1.5);

sgtitle(sprintf('Flutter Analysis for [PZT/0/90/90/0/PZT] Laminate (a=%.1f mm, b=%.1f mm)', ...
        a*1000, b*1000), 'FontSize', 16);

%% Comparison plot
figure('Position', [50, 50, 1200, 500]);

subplot(1,2,1);
plot(lambda_values, freq_Hz_damp(2,:)/1000, 'r-', 'LineWidth', 2); hold on;
plot(lambda_values, freq_Hz_no_damp(2,:)/1000, 'b--', 'LineWidth', 2);
xlabel('\lambda', 'FontSize', 12);
ylabel('Mode 2 Frequency (kHz)', 'FontSize', 12);
title('Mode 2 Frequency: With vs Without Damping', 'FontSize', 14);
legend('With Damping', 'Without Damping', 'Location', 'best');
grid on;
xline(lambda_cr, 'k--', 'LineWidth', 1.5);

subplot(1,2,2);
plot(lambda_values, damping_ratios_damp(2,:), 'r-', 'LineWidth', 2); hold on;
plot(lambda_values, damping_ratios_no_damp(2,:), 'b--', 'LineWidth', 2);
xlabel('\lambda', 'FontSize', 12);
ylabel('Mode 2 Damping Ratio', 'FontSize', 12);
title('Mode 2 Damping: With vs Without Damping', 'FontSize', 14);
legend('With Damping', 'Without Damping', 'Location', 'best');
grid on;
yline(0, 'k--', 'LineWidth', 1);
xline(lambda_cr, 'k--', 'LineWidth', 1.5);

sgtitle('Effect of Aerodynamic Damping on Flutter Characteristics', 'FontSize', 16);

%% Display results
fprintf('\n=== Results Summary ===\n');
fprintf('Geometric properties:\n');
fprintf('  Length (a): %.1f mm\n', a*1000);
fprintf('  Width (b): %.1f mm\n', b*1000);
fprintf('  Aspect ratio: %.3f\n', alpha);
fprintf('  Total thickness: %.2f mm\n', h_total*1000);
fprintf('\nMaterial properties:\n');
fprintf('  D11: %.2f N·m\n', D11);
fprintf('  D22: %.2f N·m\n', D22);
fprintf('  D12: %.2f N·m\n', D12);
fprintf('  D66: %.2f N·m\n', D66);
fprintf('  I0: %.4f kg/m²\n', I0_mass);
fprintf('\nFlutter characteristics:\n');
fprintf('  Critical λ: %.1f\n', lambda_cr);
fprintf('  Flutter frequency (with damping): %.1f Hz\n', ...
        freq_Hz_damp(2, find(lambda_values >= lambda_cr, 1)));
fprintf('  Flutter frequency (without damping): %.1f Hz\n', ...
        freq_Hz_no_damp(2, find(lambda_values >= lambda_cr, 1)));

%% Function to calculate laminate properties
function [A, B, D, As, I0, I2] = calculate_laminate_properties(layers, a, b)
    A = zeros(3); B = zeros(3); D = zeros(3); 
    As = zeros(2); I0 = 0; I2 = 0;
    
    h_tot = sum(cellfun(@(x) x.thickness, layers));
    z_cur = -h_tot/2;
    
    for i = 1:length(layers)
        L = layers{i};
        z_next = z_cur + L.thickness;
        
        if strcmp(L.type, 'comp')
            % Orthotropic composite layer
            nu21 = L.nu12 * L.E2 / L.E1;
            denom = 1 - L.nu12 * nu21;
            Q11 = L.E1 / denom;
            Q22 = L.E2 / denom;
            Q12 = L.nu12 * L.E2 / denom;
            Q66 = L.G12;
            
            Q = [Q11, Q12, 0; Q12, Q22, 0; 0, 0, Q66];
            
            % Transform for orientation angle
            theta = L.angle * pi/180;
            m = cos(theta);
            n = sin(theta);
            
            T = [m^2, n^2, 2*m*n;
                 n^2, m^2, -2*m*n;
                 -m*n, m*n, m^2-n^2];
            
            Q_trans = T \ Q * inv(T);
            
            % Shear correction factor (for FSDT)
            Qs = [L.G12, 0; 0, L.G12 * 0.6];
            Qs_trans = Qs;  % For 0/90 layers, shear transformation is identity
            
        else  % Piezoelectric layer (isotropic)
            E = L.E;
            nu = L.nu;
            denom = 1 - nu^2;
            Q = (E / denom) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
            Q_trans = Q;
            
            Qs = (E / (2*(1+nu))) * eye(2);
            Qs_trans = Qs;
        end
        
        % Update stiffness matrices
        A = A + Q_trans * (z_next - z_cur);
        B = B + Q_trans * (z_next^2 - z_cur^2) / 2;
        D = D + Q_trans * (z_next^3 - z_cur^3) / 3;
        As = As + Qs_trans * (z_next - z_cur);
        
        % Update inertia terms
        I0 = I0 + L.rho * (z_next - z_cur);
        I2 = I2 + L.rho * (z_next^3 - z_cur^3) / 3;
        
        z_cur = z_next;
    end
    
    % Apply shear correction factor (5/6 for rectangular cross-section)
    As = As * 5/6;
end
