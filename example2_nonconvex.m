%% Example 2: structured composite nonconvex coupled constraint
% Outputs: results/example2_nonconvex beside this script.

clearvars; close all; clc;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
out_dir = fullfile(script_dir, 'results', 'example2_nonconvex');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

t_final = 10;

m = shared_model();
p = proof_parameters();
[y0, init, cert] = initialization_and_certificate(m);
audit = theorem_audit(m, p, cert);

fprintf('EXAMPLE 2: STRUCTURED COMPOSITE NONCONVEX CONSTRAINT\n');
fprintf('Graph: fixed 4-regular graph with 10 nodes\n');
fprintf('lambda_2(L)=%.12f, ||L||_2=%.12f\n', m.lambda2, m.normL);
fprintf('Target multiplier lambda^*=%.4f in Lambda=[0, %.4f]\n', m.lambda_star, m.lambda_bar);
fprintf('Integration horizon T=%.1f s.\n', t_final);
fprintf('Initial action ranges: x_1 [%+.4f,%+.4f], x_2 [%+.4f,%+.4f]\n', ...
    min(init.x0(:,1)), max(init.x0(:,1)), min(init.x0(:,2)), max(init.x0(:,2)));
fprintf('Initial multiplier range: [%.4f, %.4f]\n', min(init.lambda0), max(init.lambda0));
fprintf('Initial canonical ranges: eta_1 [%.4f, %.4f], eta_2 [%+.4f, %+.4f]\n', ...
    min(init.eta0(:,1)), max(init.eta0(:,1)), min(init.eta0(:,2)), max(init.eta0(:,2)));
off_x1 = []; off_x2 = [];
for i = 1:m.N
    for j = 1:m.N
        if i ~= j
            rows = (j-1)*m.n + (1:m.n);
            off_x1(end+1) = init.X0(rows(1), i); %#ok<AGROW>
            off_x2(end+1) = init.X0(rows(2), i); %#ok<AGROW>
        end
    end
end
fprintf('Initial opponent estimate ranges: x^i_{j,1} [%+.5f, %+.5f], x^i_{j,2} [%+.5f, %+.5f]\n', ...
    min(off_x1), max(off_x1), min(off_x2), max(off_x2));
fprintf('Initial g_Sigma(x(0))=%.6e, tilde_g_Sigma(x(0),eta(0))=%.6e\n', ...
    aggregate_original_g(init.x0, m), aggregate_canonical_g(init.x0, init.eta0, m));
fprintf('Target ||omega^*||=%.6f, 1^T omega^*=%.3e\n', ...
    norm(m.omegastar), sum(m.omegastar));
fprintf('W_2(0)=%.12e, r_2=%.12e, 2r_2=%.12e\n', ...
    cert.W0, cert.r2, cert.level);
fprintf(['Certified on K_{2,2r_2}: l_theta=%.9f <= %.2f, ' ...
    'l_varphi=%.9f <= %.2f, l_phi=%.9f <= %.2f\n'], ...
    cert.ltheta_exact, m.ltheta_used, cert.lvarphi_exact, m.lvarphi_used, ...
    cert.lphi_exact, m.lphi_used);
fprintf('Uniform bounds: L_tilde_g_x=%.9g, L_G_x=%.9g, bar_rho_x=%.9g, bar_rho=%.9g, bar_Phi_star=%.9g\n', ...
    m.Lx, m.LGx, m.bar_rho_x, m.bar_rho, m.bar_Phi_star);
fprintf('Regional bounds on K_{2,2r_2}: R_x_2r=%.9g, R_eta_2r=%.9g, L_tilde_g_eta_2r=%.9g\n', ...
    cert.Rx_2r, cert.Reta_2r, cert.Leta_2r);
fprintf('Original-constraint Hessian witness at u=(0,pi/(2k)): min eigenvalue %.6e < 0\n', ...
    m.nonconvex_hessian_min);
fprintf('Target KKT residual = %.3e, canonical consistency = %.3e, Slater margin = %.6f\n', ...
    m.target_kkt_residual, m.target_consistency_residual, m.slater_margin);
disp(audit(:, {'Condition','Value','Lower','Upper','Margin','Pass'}));
writetable(audit, fullfile(out_dir, 'nonconvex_theorem_audit.csv'));

constant_name = ["q"; "E"; "a_outer"; "k_wave"; "A_amp"; "b"; "m_inner"; ...
    "lambda_star"; "lambda_bar"; "mu_G"; "l_Phi"; "l_G"; "iota"; ...
    "L_tilde_g_x"; "L_G_x"; "bar_rho_x"; "bar_rho"; "bar_Phi_star"; ...
    "R_x_2r"; "R_eta_2r"; "L_tilde_g_eta_2r"; ...
    "negative_Hessian_eigenvalue"; "slater_margin"; ...
    "W2_initial"; "r2"; "two_r2"; ...
    "certified_l_theta_K2_2r"; "certified_l_varphi_K2_2r"; "certified_l_phi_K2_2r"];
constant_value = [m.q; m.E; m.a_outer; m.wave; m.rho2_amp; m.b; m.m_inner; ...
    m.lambda_star; m.lambda_bar; m.mu_G; m.lPhi; m.lG; m.iota; ...
    m.Lx; m.LGx; m.bar_rho_x; m.bar_rho; m.bar_Phi_star; ...
    cert.Rx_2r; cert.Reta_2r; cert.Leta_2r; ...
    m.nonconvex_hessian_min; m.slater_margin; ...
    cert.W0; cert.r2; cert.level; ...
    cert.ltheta_exact; cert.lvarphi_exact; cert.lphi_exact];
model_constants = table(constant_name, constant_value, ...
    'VariableNames', {'Constant','Value'});
writetable(model_constants, fullfile(out_dir, 'nonconvex_model_constants.csv'));

agent = (1:m.N).';
equilibrium_table = table(agent, m.xstar(:,1), m.xstar(:,2), m.r(:,1), m.r(:,2), ...
    m.local_offsets, m.lambdastar, m.omegastar, ...
    'VariableNames', {'Agent','xStar1','xStar2','Reference1','Reference2', ...
    'LocalConstraintAtTarget','lambdaStar','omegaStar'});
writetable(equilibrium_table, fullfile(out_dir, 'nonconvex_equilibrium_targets.csv'));
assert(all(audit.Pass), 'Theorem-inequality check failed.');

% Export round-trip double-precision values for the Table III certificate.
fid = fopen(fullfile(out_dir, 'table_III_certificate.txt'), 'w');
assert(fid ~= -1, 'Cannot open certificate report.');
file_cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Table III numerical certificate\n');
fprintf(fid, 'MATLAB version: %s\n', version);
fprintf(fid, 'gamma = %.17g\n', p.gamma);
fprintf(fid, 'iota  = %.17g\n\n', m.iota);
for j = 1:height(audit)
    fprintf(fid, '%s\n', char(audit.Condition(j)));
    fprintf(fid, '  value = %.17g\n', audit.Value(j));
    fprintf(fid, '  lower = %.17g\n', audit.Lower(j));
    fprintf(fid, '  upper = %.17g\n', audit.Upper(j));
    fprintf(fid, '  margin = %.17g\n', audit.Margin(j));
    fprintf(fid, '  pass  = %d\n\n', audit.Pass(j));
end
clear file_cleanup

%% Canonical mirror dynamics (C-DBPD)
t_eval = unique([0, linspace(0, 0.002, 400), linspace(0.002, 0.05, 300), ...
    linspace(0.05, 0.5, 400), linspace(0.5, t_final, 1400)]);

opts = odeset('RelTol', 2e-9, 'AbsTol', 2e-11, 'MaxStep', 0.01);
[tp, yp] = ode15s(@(~, yy) proposed_rhs(yy, m, p), t_eval, y0, opts);
proposed = decode_proposed(tp, yp, m);
pm = proposed_metrics(proposed, m);

% Solver sensitivity: smaller tolerances and maximum step size.
opts_tight = odeset('RelTol', 2e-10, 'AbsTol', 2e-12, 'MaxStep', 0.005);
[tp_tight, yp_tight] = ode15s(@(~, yy) proposed_rhs(yy, m, p), t_eval, y0, opts_tight);
proposed_tight = decode_proposed(tp_tight, yp_tight, m);
pm_tight = proposed_metrics(proposed_tight, m);

post = post_audit(proposed, pm, m, cert);
disp(post);
writetable(post, fullfile(out_dir, 'nonconvex_post_simulation_audit.csv'));
assert(all(post.Pass), 'A post-simulation numerical check failed.');

%% Compare baseline and tightened solver settings
post_tight = post_audit(proposed_tight, pm_tight, m, cert);
solver_sensitivity = table(post.Check, post.Value, post_tight.Value, ...
    abs(post_tight.Value - post.Value), post.Tolerance, post.Pass, post_tight.Pass, ...
    'VariableNames', {'Check','Baseline','Tightened','AbsDifference', ...
    'Tolerance','BaselinePass','TightenedPass'});
disp(solver_sensitivity);
writetable(solver_sensitivity, fullfile(out_dir, 'nonconvex_solver_sensitivity.csv'));
assert(all(post_tight.Pass), 'A tightened-solver post-simulation numerical check failed.');

x_end = squeeze(proposed.x(end,:,:));
x_end_tight = squeeze(proposed_tight.x(end,:,:));
eta_end = squeeze(proposed.eta(end,:,:));
eta_end_tight = squeeze(proposed_tight.eta(end,:,:));
terminal_quantity = ["x terminal Frobenius difference"; ...
    "lambda terminal 2-norm difference"; ...
    "omega terminal 2-norm difference"; ...
    "eta terminal Frobenius difference"; ...
    "original-constraint terminal difference"; ...
    "canonical-constraint terminal difference"; ...
    "consistency-residual terminal difference"; ...
    "KKT-residual terminal difference"];
terminal_difference = [norm(x_end_tight - x_end, 'fro'); ...
    norm(proposed_tight.lambda(end,:) - proposed.lambda(end,:)); ...
    norm(proposed_tight.omega(end,:) - proposed.omega(end,:)); ...
    norm(eta_end_tight - eta_end, 'fro'); ...
    abs(proposed_tight.original_constraint(end) - proposed.original_constraint(end)); ...
    abs(proposed_tight.canonical_constraint(end) - proposed.canonical_constraint(end)); ...
    abs(proposed_tight.canonical_consistency(end) - proposed.canonical_consistency(end)); ...
    abs(pm_tight.kkt_terminal - pm.kkt_terminal)];
terminal_sensitivity = table(terminal_quantity, terminal_difference, ...
    'VariableNames', {'Quantity','AbsDifference'});
disp(terminal_sensitivity);
writetable(terminal_sensitivity, fullfile(out_dir, 'nonconvex_terminal_sensitivity.csv'));

make_figures(proposed, pm, m, cert, out_dir);
save(fullfile(out_dir, 'nonconvex_run_data.mat'), ...
    'm', 'p', 'init', 'cert', 'audit', 'model_constants', 'equilibrium_table', ...
    'proposed', 'pm', 'post_tight', 'solver_sensitivity', 'terminal_sensitivity', ...
    'post', '-v7.3');
fprintf('Completed. Results: %s\n', out_dir);

%% Local functions
function m = shared_model()
    m.N = 10; m.n = 2; m.d = 1;
    m.A = [ ...
        0 0 0 1 1 0 1 0 0 1; ...
        0 0 1 0 0 0 1 0 1 1; ...
        0 1 0 0 1 0 1 1 0 0; ...
        1 0 0 0 0 1 1 1 0 0; ...
        1 0 1 0 0 1 0 0 1 0; ...
        0 0 0 1 1 0 0 1 1 0; ...
        1 1 1 1 0 0 0 0 0 0; ...
        0 0 1 1 0 1 0 0 0 1; ...
        0 1 0 0 1 1 0 0 0 1; ...
        1 1 0 0 0 0 0 1 1 0];
    assert(isequal(m.A, m.A.') && all(sum(m.A, 2) == 4));
    m.L = diag(sum(m.A, 2)) - m.A;
    ev = sort(real(eig(m.L))); m.lambda2 = ev(2); m.normL = ev(end);

    % Game parameters and action bounds (moderate, low-stiffness scale)
    m.q = 400.0;
    m.beta = 0.2;
    m.h = 45.0;
    m.lambda_star = 12.0;
    m.lambda_bar = 24.0;
    m.xlo = [-1.52; -8.0];
    m.xhi = [ 1.52;  8.0];

    theta = 2*pi*(0:m.N-1)'/m.N;
    x1_target = linspace(-0.65, 0.65, m.N).';
    x2_target = 1.35*cos(theta) - 0.35*sin(2*theta);
    m.xstar = [x1_target, x2_target];

    % Local constraint offsets summing to zero
    m.local_offsets = 0.005*(-4.5:1:4.5).';
    assert(abs(sum(m.local_offsets)) < 1e-14);
    m.omegastar = pinv(m.L)*m.local_offsets;
    m.omegastar = m.omegastar - mean(m.omegastar);
    m.Romega = [eye(m.N-1); -ones(1, m.N-1)];
    m.Romega_pinv = (m.Romega.'*m.Romega)\m.Romega.';
    assert(norm(m.L*m.omegastar - m.local_offsets) < 1e-10);
    assert(abs(sum(m.omegastar)) < 1e-12);
    m.lambdastar = m.lambda_star*ones(m.N, 1);

    % Composite nonconvex constraint:
    % u_i = x_i - s_i, rho_i(x_i) = [b*u_1 + (m/2)||u||^2; A*sin(k*u_2)],
    % Phi_i(y) = E*y_1 + (a/2)||y||^2 + c_i.
    m.E = 500.0;
    m.a_outer = 6.0;
    m.mu_G = 0.2;
    m.wave = 1.5;
    m.b = m.h / m.E;              % 0.09 exact
    m.rho2_amp = 0.35;
    m.rho2_derivative_amp = m.rho2_amp * m.wave;     % 0.525
    m.rho2_curvature = m.rho2_amp * m.wave^2;        % 0.7875
    m.m_inner = 1.0e-3;
    m.eta_half = [10.0; 0.30];
    m.etastar = repmat([m.E, 0], m.N, 1);
    m.etalo = [m.E - m.eta_half(1); -m.eta_half(2)];
    m.etahi = [m.E + m.eta_half(1);  m.eta_half(2)];

    m.mu_theta = 1.0; m.mu_varphi = 1.0;
    m.mu_eta = [1.0; 1.0];
    m.mu_canonical = min(m.mu_eta);
    m.ltheta_used = 1.30;
    m.lvarphi_used = 1.05;
    m.lphi_used = 1.05;

    m.ell = m.q + m.beta*m.normL;
    m.tilde_ell = sqrt((m.q + m.beta*4)^2 + 4*m.beta^2);
    m.mu_x = m.q;
    m.r = m.xstar + (m.beta*(m.L*m.xstar) + ...
        m.lambda_star*repmat([m.h, 0], m.N, 1))/m.q;

    % Whole-domain uniform bounds on Omega x Delta
    umax = zeros(m.N, m.n);
    for i = 1:m.N
        umax(i,:) = max(abs(m.xlo.' - m.xstar(i,:)), ...
                        abs(m.xhi.' - m.xstar(i,:)));
    end
    m.umax = max(umax, [], 1);
    rho1abs = m.b*m.umax(1) + 0.5*m.m_inner*sum(m.umax.^2);
    rho_local = hypot(rho1abs, m.rho2_amp);
    m.bar_rho = sqrt(m.N)*rho_local*(1 + 1e-8);
    m.bar_Phi_star = sqrt(m.N)*hypot( ...
        m.eta_half(1)/m.a_outer, m.eta_half(2)/m.a_outer);
    base_grad = hypot(m.b + m.m_inner*m.umax(1), m.m_inner*m.umax(2));
    m.bar_rho_x = sqrt(base_grad^2 + m.rho2_derivative_amp^2)*(1 + 1e-8);
    % G(x, eta^circ) has eta_i^circ = [E; 0], so ||nabla_x G_i(x_i, eta_i^circ)|| <= E * base_grad.
    m.Lx = m.E * base_grad * (1 + 1e-8);
    % ||nabla_x(eta_i^T nabla rho_i(x_i))|| <= |eta_{i,1}|*m + |eta_{i,2}|*A*k^2.
    m.LGx = ((m.E + m.eta_half(1))*m.m_inner + ...
        m.eta_half(2)*m.rho2_curvature)*(1 + 1e-8);
    m.lG = m.h;
    m.iota = m.h^2;
    m.lPhi = m.a_outer;
    m.muPhi = m.a_outer;

    % Uniform canonical Slater witness: xhat_i = s_i - [0.05; 0]
    worst = zeros(m.N, 1);
    for i = 1:m.N
        xhat = m.xstar(i,:).' - [0.05; 0];
        assert(all(xhat > m.xlo & xhat < m.xhi));
        [rho, ~] = rho_map(xhat, i, m);
        eta1 = max(m.etalo(1), min(m.etahi(1), m.E + m.a_outer*rho(1)));
        worst(i) = eta1*rho(1) - (eta1 - m.E)^2/(2*m.a_outer) + m.local_offsets(i);
    end
    m.slater_margin = -sum(worst);
    assert(m.slater_margin > 20, 'Uniform canonical Slater witness failed.');

    m.target_interior_margin = min([ ...
        min(m.xstar - m.xlo.', [], 'all'), ...
        min(m.xhi.' - m.xstar, [], 'all'), ...
        m.lambda_star, m.lambda_bar - m.lambda_star, ...
        m.eta_half(1), m.eta_half(2)]);
    target_kkt = m.q*(m.xstar - m.r) + m.beta*(m.L*m.xstar) + ...
        m.lambda_star*repmat([m.h, 0], m.N, 1);
    m.target_kkt_residual = norm(target_kkt, 'fro');
    m.target_consistency_residual = 0.0;
    m.target_kkt_margin = 1e-10 - m.target_kkt_residual;

    % Nonconvexity witness at x_i = s_i + [0; pi/(2k)] in int(Omega_i)
    hmin = inf;
    for i = 1:m.N
        xnc = m.xstar(i,:).' + [0; pi/(2*m.wave)];
        assert(all(xnc > m.xlo & xnc < m.xhi));
        hmin = min(hmin, min(eig(original_constraint_hessian(xnc, i, m))));
    end
    m.nonconvex_hessian_min = hmin;
    assert(m.nonconvex_hessian_min < -1.0);
end

function p = proof_parameters()
    % Auxiliary parameters for Theorem 2(ii) and Table III
    p.epsilon5 = 6.8e-3;
    p.epsilon4 = 0.28;
    p.sigma    = 0.08;
    p.epsilon6 = 3.2;
    p.komega   = 6.1e3;
    p.delta3   = 3.6e-6;
    p.delta4   = p.komega * p.delta3^2;
    p.gamma    = 2200.0;
end

function [y0, init, cert] = initialization_and_certificate(m)
    N = m.N; n = m.n; k = (1:N).';
    b1 = sin(k*sqrt(2)) + 0.25*cos(k*sqrt(7));
    b1 = b1 - mean(b1); b1 = b1 / max(abs(b1));
    b2 = cos(k*sqrt(3)) - 0.2*sin(k*sqrt(5));
    b2 = b2 - mean(b2); b2 = b2 / max(abs(b2));
    [ii, jj] = ndgrid(1:N, 1:N);
    e1 = sin(ii*sqrt(3) + jj*sqrt(5)); e1 = e1 / max(abs(e1(:)));
    e2 = cos(ii*sqrt(5) - jj*sqrt(2)); e2 = e2 / max(abs(e2(:)));

    % Nonzero mean shift (+1.5e-3 in x_1) so g_Sigma(x(0)) > 0 (initially infeasible)
    x0 = m.xstar + [2.5e-3*b1 + 1.5e-3, 4.5e-3*b2];
    assert(all(x0 > m.xlo.' & x0 < m.xhi.', 'all'));
    X = zeros(n*N, N);
    for i = 1:N
        for j = 1:N
            rows = (j-1)*n + (1:n);
            if i == j
                X(rows, i) = x0(j,:).';
            else
                X(rows, i) = m.xstar(j,:).' + [1.0e-3*e1(j,i); 1.5e-3*e2(j,i)];
            end
        end
    end
    lambda0 = m.lambda_star + 5.0e-3*b1;
    omega0 = zeros(N, 1);
    zeta0 = omega0(1:N-1);
    eta0 = [m.E + 5.0e-3*b1, 3.0e-3*b2];
    assert(all(eta0(:,1) > m.etalo(1) & eta0(:,1) < m.etahi(1)));
    assert(all(eta0(:,2) > m.etalo(2) & eta0(:,2) < m.etahi(2)));

    Q = X;
    for i = 1:N
        rows = (i-1)*n + (1:n);
        Q(rows, i) = fermi_inverse(x0(i,:).', m.xlo, m.xhi, m.mu_theta);
    end
    nu0 = fermi_inverse(lambda0, 0, m.lambda_bar, m.mu_varphi);
    alpha0 = zeros(2, N);
    for i = 1:N
        alpha0(:, i) = fermi_inverse(eta0(i,:).', m.etalo, m.etahi, m.mu_eta);
    end
    y0 = [Q(:); nu0; zeta0; alpha0(:)];

    Waction = 0;
    for i = 1:N
        Waction = Waction + sum(fermi_D(m.xstar(i,:).', x0(i,:).', ...
            m.xlo, m.xhi, m.mu_theta));
    end
    Woff = 0;
    for i = 1:N
        for j = 1:N
            if i ~= j
                rows = (j-1)*n + (1:n);
                Woff = Woff + 0.5*sum((X(rows,i) - m.xstar(j,:).').^2);
            end
        end
    end
    Wlambda = sum(fermi_D(m.lambda_star*ones(N,1), lambda0, ...
        zeros(N,1), m.lambda_bar*ones(N,1), m.mu_varphi));
    Womega = 0.5*sum((omega0 - m.omegastar).^2);
    Weta = 0;
    for i = 1:N
        Weta = Weta + sum(fermi_D(m.etastar(i,:).', eta0(i,:).', ...
            m.etalo, m.etahi, m.mu_eta));
    end
    cert.W0 = Waction + Woff + Wlambda + Womega + Weta;
    cert.r2 = 1.05 * cert.W0;
    cert.level = 2.0 * cert.r2;

    % Conservative regional displacement and sensitivity bounds on K_{2,2r_2}
    cert.Rx_2r = 2.0 * sqrt(cert.r2 / m.mu_theta);
    cert.Reta_2r = 2.0 * sqrt(cert.r2 / m.mu_canonical);
    cert.Leta_2r = min(m.bar_rho + m.bar_Phi_star, ...
        m.bar_rho_x * cert.Rx_2r + cert.Reta_2r / m.muPhi);

    curv = zeros(N, n);
    for i = 1:N
        for ell = 1:n
            curv(i, ell) = fermi_level_curvature(m.xstar(i, ell), ...
                m.xlo(ell), m.xhi(ell), m.mu_theta, cert.level);
        end
    end
    cert.ltheta_exact = max(curv, [], 'all');
    cert.lvarphi_exact = fermi_level_curvature( ...
        m.lambda_star, 0, m.lambda_bar, m.mu_varphi, cert.level);
    cert.lphi1_exact = fermi_level_curvature( ...
        m.E, m.etalo(1), m.etahi(1), m.mu_eta(1), cert.level);
    cert.lphi2_exact = fermi_level_curvature( ...
        0, m.etalo(2), m.etahi(2), m.mu_eta(2), cert.level);
    cert.lphi_exact = max(cert.lphi1_exact, cert.lphi2_exact);

    assert(cert.ltheta_exact <= m.ltheta_used, ...
        'Primal mirror certificate failed: l_theta=%.9g > %.9g.', ...
        cert.ltheta_exact, m.ltheta_used);
    assert(cert.lvarphi_exact <= m.lvarphi_used, ...
        'Multiplier mirror certificate failed: l_varphi=%.9g > %.9g.', ...
        cert.lvarphi_exact, m.lvarphi_used);
    assert(cert.lphi_exact <= m.lphi_used, ...
        'Canonical mirror certificate failed: l_phi=%.9g > %.9g.', ...
        cert.lphi_exact, m.lphi_used);
    init = struct('x0', x0, 'X0', X, 'lambda0', lambda0, ...
        'omega0', omega0, 'zeta0', zeta0, 'eta0', eta0);
end

function T = theorem_audit(m, p, cert)
    % Exact implementation of Theorem 2(ii) and Appendix G
    lt = m.ltheta_used; lv = m.lvarphi_used;
    mt = m.mu_theta; mv = m.mu_varphi;
    mu2 = min([m.mu_theta, m.mu_varphi, m.mu_canonical, 1]);
    Leta_2r = cert.Leta_2r;

    lambdaSigma = m.lambda_star;
    lambdaMin   = m.lambda_star;
    mu_c  = m.mu_x + 0.5*m.mu_G*lambdaSigma;
    d_eta = lambdaMin / m.lPhi;
    c_eta = m.d * m.lambda_bar^2 * m.bar_rho_x^2;
    core_lhs = mu_c * d_eta;
    core_rhs = c_eta;

    % Sequential parameter construction (Appendix G)
    eps5_lo = 1.0 / mu_c;
    eps5_hi = d_eta / c_eta;
    mu5 = mu_c - 1.0 / p.epsilon5;
    e0  = d_eta - p.epsilon5 * c_eta;

    eps4_min = 5.0 * m.N / (2.0 * mu5);
    mu4 = mu5 - 5.0 * m.N / (2.0 * p.epsilon4);

    sigma_max = min([1.0, p.epsilon4, ...
        mu4 / (0.5*m.mu_G*lambdaSigma + 1.0/p.epsilon5)]);
    a0 = mu4 - p.sigma * (0.5*m.mu_G*lambdaSigma + 1.0/p.epsilon5);

    eps6_min = 2.0 / e0;
    e1 = e0 - 2.0 / p.epsilon6;

    bomega = p.epsilon4 * m.lG^2 * m.normL^2 / (2.0 * mv^2);
    komega_min = bomega / m.lambda2;
    domega = p.komega * m.lambda2 - bomega;
    qomega = p.epsilon4 * m.Lx^2 + p.epsilon6 * Leta_2r^2;

    % Eq. (30): delta3 compatibility & varrho2 < 1/3
    d3_cond1 = a0 / m.N - p.delta3 * m.lG * m.Lx / mv;
    d3_cond2 = e1 - p.epsilon4 * p.delta3^2 * m.lG^2 * Leta_2r^2 / (2.0 * mv^2);
    d3_cond3 = domega - 0.5 * p.komega^2 * qomega * p.delta3^2;
    varrho2  = (p.delta3 * m.lG + p.komega * p.delta3^2 * lv) / mu2;
    delta4_error = abs(p.delta4 - p.komega * p.delta3^2);

    % Gain thresholds gamma_0^{nc} and gammabar(delta3)
    tilde_beta_x = p.delta3 * m.lG * m.Lx / mv + 2.5 / p.epsilon4;
    Atilde = a0 / m.N - p.delta3 * m.lG * m.Lx / mv;

    gamma0_nc = (m.tilde_ell + (1.0/p.sigma - 1.0)*0.5*m.mu_G*lambdaSigma + ...
        (1.0 + 1.0/p.sigma)/p.epsilon5 + 2.5/p.epsilon4 + ...
        (m.ell + m.tilde_ell)^2 / (4.0 * a0)) / ...
        ((1.0 - p.sigma/p.epsilon4) * m.lambda2);

    gammabar = (m.tilde_ell + (1.0/p.sigma - 1.0)*0.5*m.mu_G*lambdaSigma + ...
        (1.0 + 1.0/p.sigma)/p.epsilon5 + tilde_beta_x + ...
        (m.ell + m.tilde_ell)^2 / (4.0 * m.N * Atilde)) / ...
        ((1.0 - p.sigma/p.epsilon4) * m.lambda2);

    % Eq. (31): iota condition evaluated at chosen gamma > gammabar
    QF = @(gam) (m.lG^2 / mt^2) * ( ...
        p.epsilon4 * (m.ell^2 / m.N + m.tilde_ell^2 + ...
        gam * m.normL^2 / (4.0 * p.sigma * m.lambda2) + ...
        0.5 * m.d * m.lambda_bar^2 * m.LGx^2) + ...
        0.5 * p.epsilon6 * c_eta );
    Qlambda = @(gam) QF(gam) + p.komega * lv * m.normL;
    iotareq = lt * (p.epsilon6 * Leta_2r^2 / p.delta3 + ...
        p.delta3 * Qlambda(p.gamma));

    % Exact Lyapunov matrix M_tilde entries
    Btilde = (1.0 - p.sigma/p.epsilon4) * p.gamma * m.lambda2 - m.tilde_ell - ...
        (1.0/p.sigma - 1.0)*0.5*m.mu_G*lambdaSigma - ...
        (1.0 + 1.0/p.sigma)/p.epsilon5 - tilde_beta_x;
    Ctilde = p.delta3 * m.iota / lt - p.epsilon6 * Leta_2r^2 - ...
        p.delta3^2 * Qlambda(p.gamma);
    Dtilde = p.delta4 * m.lambda2 - ...
        p.epsilon4 * p.delta3^2 * m.lG^2 * m.normL^2 / (2.0 * mv^2) - ...
        0.5 * p.delta4^2 * (p.epsilon4 * m.Lx^2 + p.epsilon6 * Leta_2r^2);
    Etilde = d_eta - p.epsilon5 * c_eta - ...
        p.epsilon4 * p.delta3^2 * m.lG^2 * Leta_2r^2 / (2.0 * mv^2) - ...
        2.0 / p.epsilon6;
    Schur = Btilde - (m.ell + m.tilde_ell)^2 / (4.0 * m.N * Atilde);

    cross = -(m.ell + m.tilde_ell) / (2.0 * sqrt(m.N));
    M = zeros(5);
    M(1:2, 1:2) = [Atilde, cross; cross, Btilde];
    M(3,3) = Ctilde; M(4,4) = Dtilde; M(5,5) = Etilde;
    minM = min(eig(M));

    canonical_curvature = (m.E - m.eta_half(1))*m.m_inner - ...
        m.eta_half(2)*m.rho2_curvature;

    cond = { ...
        'graph algebraic connectivity'; ...
        'uniform canonical Slater margin'; ...
        'target interior margin'; ...
        'target KKT construction margin'; ...
        'Delta subset int(E_plus): canonical curvature'; ...
        'negative original-constraint Hessian eigenvalue'; ...
        'Bregman sublevel r2 > W2(0)'; ...
        'regional l_theta <= 1.30'; ...
        'regional l_varphi <= 1.05'; ...
        'regional l_phi <= 1.05'; ...
        'core small gain (12): mu_c*d_eta > c_eta'; ...
        'epsilon5 in (1/mu_c, d_eta/c_eta)'; ...
        'epsilon4 > 5N/(2mu5)'; ...
        'sigma in (0, min{1,eps4,mu4/(0.5*mu_G*lam_Sigma+1/eps5)})'; ...
        'epsilon6 > 2/e0'; ...
        'komega > bomega/lambda2'; ...
        'delta3 condition 1 (Eq. 30)'; ...
        'delta3 condition 2 (Eq. 30)'; ...
        'delta3 condition 3 (Eq. 30)'; ...
        'delta4 = komega*delta3^2'; ...
        'varrho2 < 1/3 (Eq. 30)'; ...
        'gamma > gammabar(delta3)'; ...
        'iota > Eq. (31) threshold'; ...
        'A_tilde > 0'; ...
        'Schur margin > 0'; ...
        'C_tilde/delta3 > 0'; ...
        'D_tilde/delta3^2 > 0'; ...
        'E_tilde > 0'; ...
        'lambda_min(M_tilde) > 0'};

    val = [m.lambda2; m.slater_margin; m.target_interior_margin; ...
        m.target_kkt_margin; canonical_curvature; ...
        -m.nonconvex_hessian_min; cert.r2; ...
        cert.ltheta_exact; cert.lvarphi_exact; cert.lphi_exact; ...
        core_lhs; p.epsilon5; p.epsilon4; p.sigma; p.epsilon6; ...
        p.komega; d3_cond1; d3_cond2; d3_cond3; delta4_error; ...
        varrho2; p.gamma; m.iota; Atilde; Schur; ...
        Ctilde/p.delta3; Dtilde/p.delta3^2; Etilde; minM];

    lo = [0; 0; 0; 0; m.mu_G; 0; cert.W0; ...
        0; 0; 0; core_rhs; eps5_lo; eps4_min; 0; eps6_min; komega_min; ...
        0; 0; 0; -inf; 0; gammabar; iotareq; 0; 0; 0; 0; 0; 0];
    hi = [inf; inf; inf; inf; inf; inf; inf; ...
        m.ltheta_used; m.lvarphi_used; m.lphi_used; inf; eps5_hi; inf; sigma_max; ...
        inf; inf; inf; inf; inf; 1e-24; 1/3; inf; inf; inf; inf; inf; inf; inf; inf];

    pass = val > lo & val < hi;
    pass(20) = delta4_error <= 1e-24;

    margin = min((val - lo)./max([abs(val), abs(lo), ones(size(val))*1e-14], [], 2), ...
                 (hi - val)./max([abs(val), abs(hi), ones(size(val))*1e-14], [], 2));
    margin(20) = max(0, 1 - delta4_error/1e-24);
    infinite_hi = ~isfinite(hi);
    margin(infinite_hi) = (val(infinite_hi) - lo(infinite_hi))./ ...
        max([abs(val(infinite_hi)), abs(lo(infinite_hi)), ...
        ones(sum(infinite_hi), 1)*1e-14], [], 2);

    T = table(string(cond), val, lo, hi, margin, pass, ...
        'VariableNames', {'Condition','Value','Lower','Upper','Margin','Pass'});

    fprintf(['Table-III certificate: core %.12e > %.12e; eps5=%.12e in ' ...
        '(%.12e, %.12e); eps4_min=%.12e; sigma_max=%.12e; ' ...
        'eps6_min=%.12e; komega_min=%.12e\n'], ...
        core_lhs, core_rhs, p.epsilon5, eps5_lo, eps5_hi, eps4_min, sigma_max, ...
        eps6_min, komega_min);
    fprintf(['gamma0_nc=%.12e, gammabar=%.12e, gamma=%.12e; delta3=%.12e, ' ...
        'delta4=%.12e, varrho2=%.12e, iota_req=%.12e (iota=%.12e)\n'], ...
        gamma0_nc, gammabar, p.gamma, p.delta3, p.delta4, varrho2, iotareq, m.iota);
    fprintf(['Lyapunov margins: A=%.12e, Schur=%.12e, C/d3=%.12e, ' ...
        'D/d3^2=%.12e, E=%.12e, lam_min(M)=%.12e\n'], ...
        Atilde, Schur, Ctilde/p.delta3, Dtilde/p.delta3^2, Etilde, minM);
end

function dy = proposed_rhs(y, m, p)
    N = m.N; n = m.n; nq = n*N*N; nw = N - 1;
    Q = reshape(y(1:nq), n*N, N);
    nu = y(nq+(1:N));
    zeta = y(nq+N+(1:nw));
    omega = m.Romega*zeta;
    alpha_start = nq + N + nw;
    alpha = reshape(y(alpha_start+(1:2*N)), 2, N);
    [X, x] = reconstruct(Q, m);
    lambda = fermi_map(nu, 0, m.lambda_bar, m.mu_varphi);
    eta = zeros(2, N);
    for i = 1:N
        eta(:, i) = fermi_map(alpha(:, i), m.etalo, m.etahi, m.mu_eta);
    end
    dQ = zeros(size(Q)); dnu = zeros(N, 1); dalpha = zeros(2, N);
    for i = 1:N
        rowsi = (i-1)*n + (1:n); xi = x(i,:).';
        Fi = m.q*(xi - m.r(i,:).');
        for j = find(m.A(i,:))
            rowsj = (j-1)*n + (1:n);
            Fi = Fi + m.beta*(xi - X(rowsj, i));
        end
        [rho, J] = rho_map(xi, i, m);
        own_cons = zeros(n, 1);
        for j = find(m.A(i,:))
            own_cons = own_cons + (xi - X(rowsi, j));
        end
        dQ(rowsi, i) = -Fi - J.'*eta(:, i)*lambda(i) - p.gamma*own_cons;
        for player = 1:N
            if player == i, continue; end
            rows = (player-1)*n + (1:n); v = zeros(n, 1);
            for j = find(m.A(i,:))
                v = v + (X(rows, i) - X(rows, j));
            end
            dQ(rows, i) = -p.gamma*v;
        end
        gradstar = [(eta(1,i) - m.E)/m.a_outer; eta(2,i)/m.a_outer];
        phistar = ((eta(1,i) - m.E)^2 + eta(2,i)^2)/(2*m.a_outer) - m.local_offsets(i);
        dnu(i) = eta(:, i).'*rho - phistar - (m.L(i,:)*omega);
        dalpha(:, i) = lambda(i)*(rho - gradstar);
    end
    domega = m.L*lambda;
    dzeta = m.Romega_pinv*domega;
    dy = [dQ(:); dnu; dzeta; dalpha(:)];
end

function [X, x] = reconstruct(Q, m)
    X = Q; x = zeros(m.N, m.n);
    for i = 1:m.N
        rows = (i-1)*m.n + (1:m.n);
        xi = fermi_map(Q(rows, i), m.xlo, m.xhi, m.mu_theta);
        X(rows, i) = xi; x(i,:) = xi.';
    end
end

function tr = decode_proposed(t, y, m)
    nt = numel(t); N = m.N; n = m.n; nq = n*N*N; nw = N - 1;
    tr.t = t; tr.x = zeros(nt, N, n); tr.lambda = zeros(nt, N);
    tr.omega = zeros(nt, N); tr.eta = zeros(nt, N, 2);
    alpha_start = nq + N + nw;
    tr.esterr = zeros(nt, 1); tr.esterr_agent = zeros(nt, N);
    tr.actionerr = zeros(nt, 1); tr.lambdaerr = zeros(nt, 1);
    tr.omegaerr = zeros(nt, 1); tr.etaerr = zeros(nt, 1);
    tr.totalerr = zeros(nt, 1); tr.W2 = zeros(nt, 1);
    tr.original_constraint = zeros(nt, 1); tr.canonical_constraint = zeros(nt, 1);
    tr.canonical_consistency = zeros(nt, 1);
    for k = 1:nt
        Q = reshape(y(k, 1:nq), n*N, N); [X, x] = reconstruct(Q, m);
        tr.x(k,:,:) = reshape(x, 1, N, n);
        la = fermi_map(y(k, nq+(1:N)).', 0, m.lambda_bar, m.mu_varphi).';
        tr.lambda(k,:) = la;
        zeta = y(k, nq+N+(1:nw)).';
        tr.omega(k,:) = (m.Romega*zeta).';
        al = reshape(y(k, alpha_start+(1:2*N)), 2, N);
        eta = zeros(2, N);
        for i = 1:N
            eta(:, i) = fermi_map(al(:, i), m.etalo, m.etahi, m.mu_eta);
        end
        tr.eta(k,:,:) = reshape(eta.', 1, N, 2);
        go = 0; gc = 0; cc = 0;
        Wact = 0; Woff = 0; Weta = 0;
        for i = 1:N
            [rho, ~] = rho_map(x(i,:).', i, m);
            go = go + m.E*rho(1) + 0.5*m.a_outer*(rho.'*rho) + m.local_offsets(i);
            phistar = ((eta(1,i) - m.E)^2 + eta(2,i)^2)/(2*m.a_outer) - m.local_offsets(i);
            gc = gc + eta(:, i).'*rho - phistar;
            etatrue = [m.E + m.a_outer*rho(1); m.a_outer*rho(2)];
            cc = cc + norm(eta(:, i) - etatrue)^2;
            Wact = Wact + sum(fermi_D(m.xstar(i,:).', x(i,:).', m.xlo, m.xhi, m.mu_theta));
            Weta = Weta + sum(fermi_D(m.etastar(i,:).', eta(:,i), m.etalo, m.etahi, m.mu_eta));
            for j = 1:N
                if j ~= i
                    rows = (j-1)*n + (1:n);
                    Woff = Woff + 0.5*sum((X(rows, i) - m.xstar(j,:).').^2);
                end
            end
        end
        Wlam = sum(fermi_D(m.lambda_star*ones(N,1), la.', ...
            zeros(N,1), m.lambda_bar*ones(N,1), m.mu_varphi));
        Wom = 0.5*sum((tr.omega(k,:).' - m.omegastar).^2);
        tr.W2(k) = Wact + Woff + Wlam + Wom + Weta;
        tr.original_constraint(k) = go;
        tr.canonical_constraint(k) = gc;
        tr.canonical_consistency(k) = sqrt(cc);
        xstack = reshape(x.', [], 1);
        Xtruth = repmat(xstack, 1, N);
        tr.esterr(k) = norm(X - Xtruth, 'fro');
        for agent = 1:N
            tr.esterr_agent(k, agent) = norm(X(:, agent) - xstack);
        end
        tr.actionerr(k) = norm(x - m.xstar, 'fro');
        tr.lambdaerr(k) = norm(la - m.lambda_star);
        tr.omegaerr(k) = norm(tr.omega(k,:).' - m.omegastar);
        tr.etaerr(k) = norm(eta - m.etastar.', 'fro');
        tr.totalerr(k) = tr.actionerr(k) + tr.esterr(k) + tr.lambdaerr(k) + ...
            tr.omegaerr(k) + tr.etaerr(k);
    end
end

function met = proposed_metrics(tr, m)
    met.normalized_error = tr.totalerr / max(tr.totalerr(1), realmin);
    met.log_error = log(max(met.normalized_error, 1e-14));
    met.consensus = max(abs(tr.lambda - mean(tr.lambda, 2)), [], 2);
    met.omega_error = sqrt(sum((tr.omega - m.omegastar.').^2, 2));
    met.log_eta = log(max(tr.etaerr / max(tr.etaerr(1), realmin), 1e-14));
    met.log_consistency = log(max(tr.canonical_consistency / ...
        max(tr.canonical_consistency(1), realmin), 1e-14));
    idx = tr.t >= 0.5 & tr.t <= 8.0 & met.normalized_error > 1e-11;
    c = polyfit(tr.t(idx), met.log_error(idx), 1);
    met.fit_rate = -c(1); met.log_fit = polyval(c, tr.t);
    xend = squeeze(tr.x(end,:,:)); laend = tr.lambda(end,:).';
    etaend = squeeze(tr.eta(end,:,:));
    R = zeros(m.N, m.n);
    for i = 1:m.N
        [~, J] = rho_map(xend(i,:).', i, m);
        Fi = m.q*(xend(i,:).' - m.r(i,:).') + m.beta*(m.L(i,:)*xend).';
        R(i,:) = (Fi + J.'*etaend(i,:).'*laend(i)).';
    end
    met.kkt_terminal = norm(R, 'fro');
    met.max_W2 = max(tr.W2);
end

function T = post_audit(pr, pm, m, cert)
    names = ["terminal normalized total error"; ...
        "terminal estimate disagreement"; ...
        "terminal multiplier disagreement"; ...
        "terminal active original constraint"; ...
        "terminal active canonical constraint"; ...
        "terminal canonical consistency"; ...
        "terminal KKT residual"; ...
        "positive fitted exponential rate"; ...
        "zero-sum omega invariant"; ...
        "terminal omega-target residual"; ...
        "confinement max W2(t) <= 2*r2"; ...
        "negative original-constraint Hessian eigenvalue"];
    values = [pm.normalized_error(end); pr.esterr(end); pm.consensus(end); ...
        abs(pr.original_constraint(end)); abs(pr.canonical_constraint(end)); ...
        pr.canonical_consistency(end); pm.kkt_terminal; pm.fit_rate; ...
        max(abs(sum(pr.omega, 2))); pm.omega_error(end); pm.max_W2; ...
        -m.nonconvex_hessian_min];
    tolerances = [1e-9; 1e-10; 1e-9; 1e-10; 1e-10; 1e-9; 1e-7; 0; 1e-12; 1e-9; cert.level; 0];
    passes = [values(1:7) <= tolerances(1:7); values(8) > 0; ...
        values(9) <= tolerances(9); values(10) <= tolerances(10); ...
        values(11) <= tolerances(11); values(12) > 0];
    T = table(names, values, tolerances, passes, ...
        'VariableNames', {'Check','Value','Tolerance','Pass'});
end

function g = original_g(x, i, m)
    rho = rho_map(x, i, m);
    g = m.E*rho(1) + 0.5*m.a_outer*(rho.'*rho) + m.local_offsets(i);
end

function [rho, J] = rho_map(x, i, m)
    u = x - m.xstar(i,:).';
    rho = [m.b*u(1) + 0.5*m.m_inner*(u(1)^2 + u(2)^2); ...
           m.rho2_amp*sin(m.wave*u(2))];
    J = [m.b + m.m_inner*u(1), m.m_inner*u(2); ...
         0, m.rho2_derivative_amp*cos(m.wave*u(2))];
end

function H = original_constraint_hessian(x, i, m)
    [rho, J] = rho_map(x, i, m);
    eta = [m.E + m.a_outer*rho(1); m.a_outer*rho(2)];
    u = x - m.xstar(i,:).';
    H = J.'*(m.a_outer*eye(2))*J + eta(1)*m.m_inner*eye(2);
    H(2,2) = H(2,2) + eta(2)*(-m.rho2_curvature*sin(m.wave*u(2)));
end

function G = aggregate_original_g(x, m)
    G = 0;
    for i = 1:m.N
        G = G + original_g(x(i,:).', i, m);
    end
end

function G = aggregate_canonical_g(x, eta, m)
    G = 0;
    for i = 1:m.N
        rho = rho_map(x(i,:).', i, m);
        phistar = ((eta(i,1) - m.E)^2 + eta(i,2)^2)/(2*m.a_outer) - m.local_offsets(i);
        G = G + eta(i,:)*rho - phistar;
    end
end

function make_figures(pr, pm, m, cert, out_dir) %#ok<INUSD>
    cols = agent_colors();
    blue  = [0.0000, 0.0000, 1.0000];
    red   = [1.0000, 0.0000, 0.0000];
    green = [0.0000, 0.7200, 0.0000];
    nt = numel(pr.t);

    aggregate_action_error = zeros(nt, 1);
    for i = 1:m.N
        dx = squeeze(pr.x(:, i, :)) - m.xstar(i, :);
        aggregate_action_error = aggregate_action_error + sum(dx.^2, 2);
    end
    aggregate_estimate_error = sum(pr.esterr_agent.^2, 2);


    % =========================================================================
    % Fig. 6: Primal and decision-estimate squared errors (matches convex Fig. 3)
    % =========================================================================
    f = new_manuscript_figure();
    tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    t_ms = 1e3 * pr.t; % Time in s x 10^{-3}

    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, t_ms, 1e4 * aggregate_action_error, ...
        'Color', blue, 'LineStyle', '-', 'LineWidth', 5.0, ...
        'DisplayName', '$\sum_{i=1}^{10}\|x_i-x_i^\circ\|^2$');
    xlim(ax, [0, 1.0]); ylim(ax, [0, 1.3]);
    xticks(ax, [0, 0.5, 1.0]); yticks(ax, [0, 0.4, 0.8, 1.2]);
    xlabel(ax, 'Time $t$ (s) $\times 10^{-3}$', 'Interpreter', 'latex', 'FontSize', 55);
    ylabel(ax, '$\sum_{i=1}^{10}\|x_i-x_i^\circ\|^2 \times 10^{-4}$', ...
        'Interpreter', 'latex', 'FontSize', 52);
    lgd = legend(ax, 'Location', 'northeast', 'Interpreter', 'latex', ...
        'FontSize', 48, 'Box', 'off');
    lgd.ItemTokenSize = [30, 18];

    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, t_ms, 1e3 * aggregate_estimate_error, ...
        'Color', blue, 'LineStyle', '-', 'LineWidth', 5.0, ...
        'DisplayName', '$\sum_{i,j=1}^{10}\|x_j^i-x_j\|^2$');
    xlim(ax, [0, 1.0]); ylim(ax, [0, 1.35]);
    xticks(ax, [0, 0.5, 1.0]); yticks(ax, [0, 0.4, 0.8, 1.2]);
    xlabel(ax, 'Time $t$ (s) $\times 10^{-3}$', 'Interpreter', 'latex', 'FontSize', 55);
    ylabel(ax, '$\sum_{i,j=1}^{10}\|x_j^i-x_j\|^2 \times 10^{-3}$', ...
        'Interpreter', 'latex', 'FontSize', 52);
    lgd2 = legend(ax, 'Location', 'northeast', 'Interpreter', 'latex', ...
        'FontSize', 48, 'Box', 'off');
    lgd2.ItemTokenSize = [30, 18];

    format_axes(f);
    save_figure_bundle(f, out_dir, 'nonconvex_tra');
    close(f);

    % =========================================================================
    % Fig. 7: Multiplier deviations and auxiliary states (matches convex Fig. 4)
    % =========================================================================
    f = new_manuscript_figure();
    tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax = nexttile(tl); hold(ax, 'on');
    hlam = gobjects(m.N, 1);
    lambda_err = 1e3 * (pr.lambda - m.lambda_star);
    for i = 1:m.N
        hlam(i) = plot(ax, pr.t, lambda_err(:, i), 'Color', cols(i,:), ...
            'LineStyle', '-', 'LineWidth', 5.0, ...
            'DisplayName', sprintf('Node %d', i));
    end
    href = yline(ax, 0, 'Color', [0 0 0], 'LineStyle', '-.', 'LineWidth', 4.0, ...
        'DisplayName', '$\lambda_i=\tilde{\lambda}^{\circ}$');
    xlim(ax, [0, 5]); ylim(ax, [-5.2, 5.2]);
    xticks(ax, [0, 2, 4]); yticks(ax, [-5, 0, 5]);
    xlabel(ax, 'Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', 55);
    ylabel(ax, '$10^3(\lambda_i-\tilde{\lambda}^{\circ})$', ...
        'Interpreter', 'latex', 'FontSize', 55);
    lgd = legend(ax, [hlam; href], 'Interpreter', 'latex', ...
        'FontSize', 38, 'NumColumns', 1, 'Box', 'off');
    lgd.ItemTokenSize = [24, 14];
    lgd.Position = [0.338, 0.238, 0.145, 0.715];

    ax = nexttile(tl); hold(ax, 'on');
    omega_scaled = 1e3 * pr.omega;
    for i = 1:m.N
        plot(ax, pr.t, omega_scaled(:, i), 'Color', cols(i,:), ...
            'LineStyle', '-', 'LineWidth', 5.0);
    end
    plot(ax, pr.t, sum(omega_scaled, 2), 'Color', [0 0 0], ...
        'LineStyle', '-.', 'LineWidth', 4.0);
    xlim(ax, [0, 5]); ylim(ax, [-6.2, 6.2]);
    xticks(ax, [0, 2, 4]); yticks(ax, [-6, -3, 0, 3, 6]);
    xlabel(ax, 'Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', 55);
    ylabel(ax, '$10^3\omega_i$', 'Interpreter', 'latex', 'FontSize', 55);

    format_axes(f);
    save_figure_bundle(f, out_dir, 'nonconvex_aux');
    close(f);

    % =========================================================================
    % Fig. 8: Log error decay & constraint residuals (matches convex Fig. 5)
    % =========================================================================
    f = new_manuscript_figure();
    tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax = nexttile(tl); hold(ax, 'on');
    h1 = plot(ax, pr.t, pm.log_error, 'LineStyle', '-', 'LineWidth', 5.0, ...
        'Color', blue, 'DisplayName', '$e_{\Sigma,2}$');
    h2 = plot(ax, pr.t, pm.log_eta, 'LineStyle', '-', 'LineWidth', 5.0, ...
        'Color', red, 'DisplayName', '$e_{\eta}=\|\eta-\eta^{\circ}\|_{F}$');
    h3 = plot(ax, pr.t, pm.log_consistency, 'LineStyle', '--', 'LineWidth', 4.5, ...
        'Color', green, 'DisplayName', '$e_{\rm c}=\|\eta-\nabla\Phi(\rho(x))\|_{F}$');
    xlim(ax, [0, 8]); ylim(ax, [-18, 0]);
    xticks(ax, [0, 2, 4, 6, 8]); yticks(ax, [-18, -12, -6, 0]);
    xlabel(ax, 'Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', 55);
    ylabel(ax, '$\ln(e(t)/e(0))$', 'Interpreter', 'latex', 'FontSize', 55);
    lgd = legend(ax, [h1 h2 h3], 'Location', 'southwest', ...
        'Interpreter', 'latex', 'FontSize', 42, 'Box', 'off');
    lgd.ItemTokenSize = [28, 16];

    ax = nexttile(tl); hold(ax, 'on');
    h1 = plot(ax, pr.t, pr.original_constraint, 'LineStyle', '-', ...
        'LineWidth', 5.0, 'Color', blue, 'DisplayName', 'original $g_\Sigma(x)$');
    h2 = plot(ax, pr.t, pr.canonical_constraint, 'LineStyle', '--', ...
        'LineWidth', 4.5, 'Color', red, 'DisplayName', 'canonical $\tilde{g}_\Sigma(x,\eta)$');
    yline(ax, 0, 'Color', [0 0 0], 'LineStyle', '-.', 'LineWidth', 3.0, ...
        'HandleVisibility', 'off');
    xlim(ax, [0, 0.4]); ylim(ax, [-0.06, 0.72]);
    xticks(ax, [0, 0.2, 0.4]); yticks(ax, [0, 0.2, 0.4, 0.6]);
    xlabel(ax, 'Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', 55);
    ylabel(ax, 'Constraint residual', 'Interpreter', 'latex', 'FontSize', 55);
    lgd2 = legend(ax, [h1 h2], 'Location', 'northeast', ...
        'Interpreter', 'latex', 'FontSize', 44, 'Box', 'off');
    lgd2.ItemTokenSize = [28, 16];

    format_axes(f);
    save_figure_bundle(f, out_dir, 'nonconvex_exp');
    close(f);
end

function f = new_manuscript_figure()
    % Exact 2:1 widescreen canvas matching convex_aux.fig / nonconvex_com.fig
    f = figure('Color', 'w', 'Units', 'pixels', 'Position', [1 1 2048 1030], ...
        'Visible', 'off');
    set(f, 'PaperUnits', 'points', 'PaperPosition', [0 0 1485 750], ...
        'PaperSize', [1485 750]);
end

function cols = agent_colors()
    % Exact 10-node RGB palette extracted from convex_aux.fig (Fig. 4)
    cols = [ ...
        0.0000, 0.0000, 1.0000; ...
        1.0000, 0.0000, 0.0000; ...
        0.0000, 0.6275, 0.0000; ...
        1.0000, 0.5490, 0.0000; ...
        0.5020, 0.0000, 1.0000; ...
        0.0000, 0.7216, 0.8510; ...
        1.0000, 0.0000, 0.6588; ...
        0.5451, 0.2706, 0.0745; ...
        0.5020, 0.5020, 0.0000; ...
        0.2510, 0.2510, 0.2510];
end

function format_axes(f)
    ax = findall(f, 'Type', 'axes');
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 50, 'LineWidth', 0.90, ...
        'TickLabelInterpreter', 'latex', 'Box', 'on', ...
        'XGrid', 'off', 'YGrid', 'off', 'YMinorGrid', 'off');
end

function save_figure_bundle(f, out_dir, stem)
    set(f, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    savefig(f, fullfile(out_dir, [stem '.fig']));
    exportgraphics(f, fullfile(out_dir, [stem '.png']), 'Resolution', 200);
    try
        exportgraphics(f, fullfile(out_dir, [stem '.eps']), 'ContentType', 'vector');
    catch
        print(f, fullfile(out_dir, [stem '.eps']), '-depsc', '-painters');
    end
end

function x = fermi_map(z, lo, hi, mu)
    z = z(:); lo = expand(lo, numel(z)); hi = expand(hi, numel(z)); mu = expand(mu, numel(z));
    width = hi - lo; kk = mu.*width/4; a = z./kk; s = zeros(size(z)); pos = a >= 0;
    s(pos) = 1./(1 + exp(-min(a(pos), 700)));
    ea = exp(max(a(~pos), -700)); s(~pos) = ea./(1 + ea);
    x = lo + width.*s;
end

function z = fermi_inverse(x, lo, hi, mu)
    x = x(:); lo = expand(lo, numel(x)); hi = expand(hi, numel(x)); mu = expand(mu, numel(x));
    width = hi - lo; kk = mu.*width/4; p = min(max((x - lo)./width, 1e-14), 1 - 1e-14);
    z = kk.*log(p./(1 - p));
end

function D = fermi_D(ref, current, lo, hi, mu)
    ref = ref(:); current = current(:); lo = expand(lo, numel(ref)); hi = expand(hi, numel(ref));
    mu = expand(mu, numel(ref)); width = hi - lo; kk = mu.*width/4;
    pr = min(max((ref - lo)./width, 1e-14), 1 - 1e-14);
    pc = min(max((current - lo)./width, 1e-14), 1 - 1e-14);
    D = kk.*width.*(pr.*log(pr./pc) + (1 - pr).*log((1 - pr)./(1 - pc)));
end

function c = fermi_level_curvature(ref, lo, hi, mu, level)
    width = hi - lo;
    pref = min(max((ref - lo)/width, 1e-12), 1 - 1e-12);
    scale = mu*width^2/4;
    Dp = @(p) scale*(pref*log(pref./p) + (1 - pref)*log((1 - pref)./(1 - p))) - level;
    left_hi = max(pref - 1e-12, 2e-14);
    right_lo = min(pref + 1e-12, 1 - 2e-14);
    if Dp(1e-14) <= 0 || Dp(1 - 1e-14) <= 0
        error('Requested Bregman sublevel reaches too close to a boundary.');
    end
    pleft = fzero(Dp, [1e-14, left_hi]);
    pright = fzero(Dp, [right_lo, 1 - 1e-14]);
    c = max(mu/(4*pleft*(1 - pleft)), mu/(4*pright*(1 - pright)));
end

function v = expand(v, n)
    if isscalar(v), v = repmat(v, n, 1); else, v = v(:); end
end
