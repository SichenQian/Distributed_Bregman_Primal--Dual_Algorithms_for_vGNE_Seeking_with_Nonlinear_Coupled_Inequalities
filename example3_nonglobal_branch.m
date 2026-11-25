%% Example 3: distributed global vGNE recovery versus a nonglobal KKT branch
% Compares direct DBPD (5) on the original nonconvex constraint against
% proposed C-DBPD (11) from identical (x(0), X(0), lambda(0), omega(0))
% initialization in a four-agent coupled game over an undirected four-cycle.
% Outputs: results/example3_nonglobal beside this script.

clearvars; close all; clc;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
out_dir = fullfile(script_dir, 'results', 'example3_nonglobal');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

d = build_model();
init = build_initialization(d);

fprintf('EXAMPLE 3: DISTRIBUTED GLOBAL RECOVERY VS. NONGLOBAL KKT BRANCH\n');
fprintf('4-agent coupled game on 4-cycle: lambda_2(L)=%.4f, ||L||=%.4f\n', d.lambda2, d.normL);
fprintf('Global vGNE x^circ = [%.2f, %.2f, %.2f, %.2f], lambda^circ = %.2f\n', ...
    d.x_global(1), d.x_global(2), d.x_global(3), d.x_global(4), d.lam_global);
fprintf('Nonglobal KKT x_KKT = [%.2f, %.2f, %.2f, %.2f], lambda_KKT = %.2f\n', ...
    d.x_local(1), d.x_local(2), d.x_local(3), d.x_local(4), d.lam_local);
fprintf('q = [%.6f, %.6f, %.6f, %.6f]\n', d.q(1), d.q(2), d.q(3), d.q(4));
fprintf('r = [%.6f, %.6f, %.6f, %.6f]\n', d.r(1), d.r(2), d.r(3), d.r(4));
fprintf('eta^circ = [%.6f, %.6f, %.6f, %.6f] in Delta_i = [%.2f, %.2f]\n', ...
    d.eta_star(1), d.eta_star(2), d.eta_star(3), d.eta_star(4), d.etalo, d.etahi);
fprintf('Uniform canonical Slater margin at xhat=0: %.6f > 0\n', d.slater_margin);
fprintf('Global VI test at x_KKT: <F(x_KKT), x^circ - x_KKT> = %.6f < 0\n', d.vi_at_local);
fprintf('Exact unilateral global BR gap: at x_KKT = %.6f, at x^circ = %.3e\n', ...
    d.gap_at_local, d.gap_at_global);
fprintf('Exact KKT residuals: global = %.3e, local = %.3e\n', ...
    d.kkt_at_global, d.kkt_at_local);

t_final = 80.0;
t_eval = unique([0, linspace(0, 5, 400), linspace(5, t_final, 1800)]).';

opts = odeset('RelTol', 2e-10, 'AbsTol', 2e-12, 'MaxStep', 0.02);
[td, yd] = ode15s(@(~, y) rhs_direct_dbpd(y, d), t_eval, init.y0_direct, opts);
[tc, yc] = ode15s(@(~, y) rhs_cdbpd(y, d), t_eval, init.y0_cdbpd, opts);

direct = decode_direct(td, yd, d);
cdbpd  = decode_cdbpd(tc, yc, d);

% Tightened solver sensitivity check
opts_tight = odeset('RelTol', 2e-11, 'AbsTol', 2e-13, 'MaxStep', 0.01);
[td_t, yd_t] = ode15s(@(~, y) rhs_direct_dbpd(y, d), t_eval, init.y0_direct, opts_tight);
[tc_t, yc_t] = ode15s(@(~, y) rhs_cdbpd(y, d), t_eval, init.y0_cdbpd, opts_tight);
direct_t = decode_direct(td_t, yd_t, d);
cdbpd_t  = decode_cdbpd(tc_t, yc_t, d);

fprintf('\n--- TERMINAL COMPARISON (T = %.1f s, IDENTICAL INITIALIZATION) ---\n', t_final);
fprintf(['Direct DBPD (5):     x=[%+.6f %+.6f %+.6f %+.6f], ||x-x_KKT||=%.3e, ' ...
    '||lam-lam_KKT||=%.3e, KKT=%.3e, G_BR=%.6f, max|1^T om|=%.2e\n'], ...
    direct.x(end,1), direct.x(end,2), direct.x(end,3), direct.x(end,4), ...
    norm(direct.x(end,:).' - d.x_local), norm(direct.lambda(end,:).' - d.lam_local), ...
    direct.kkt(end), direct.gap(end), max(abs(sum(direct.omega, 2))));
fprintf(['Proposed C-DBPD (11): x=[%+.6f %+.6f %+.6f %+.6f], ||x-x^circ||=%.3e, ' ...
    '||lam-lam^circ||=%.3e, KKT=%.3e, e_c=%.3e, G_BR=%.3e, max|1^T om|=%.2e\n'], ...
    cdbpd.x(end,1), cdbpd.x(end,2), cdbpd.x(end,3), cdbpd.x(end,4), ...
    norm(cdbpd.x(end,:).' - d.x_global), norm(cdbpd.lambda(end,:).' - d.lam_global), ...
    cdbpd.kkt(end), cdbpd.ec(end), cdbpd.gap(end), max(abs(sum(cdbpd.omega, 2))));

sens_qty = ["Direct DBPD terminal x 2-norm difference"; ...
    "Direct DBPD terminal lambda 2-norm difference"; ...
    "Direct DBPD terminal KKT difference"; ...
    "Direct DBPD terminal G_BR difference"; ...
    "C-DBPD terminal x 2-norm difference"; ...
    "C-DBPD terminal lambda 2-norm difference"; ...
    "C-DBPD terminal eta 2-norm difference"; ...
    "C-DBPD terminal KKT difference"; ...
    "C-DBPD terminal consistency difference"; ...
    "C-DBPD terminal G_BR difference"];
sens_diff = [norm(direct_t.x(end,:) - direct.x(end,:)); ...
    norm(direct_t.lambda(end,:) - direct.lambda(end,:)); ...
    abs(direct_t.kkt(end) - direct.kkt(end)); ...
    abs(direct_t.gap(end) - direct.gap(end)); ...
    norm(cdbpd_t.x(end,:) - cdbpd.x(end,:)); ...
    norm(cdbpd_t.lambda(end,:) - cdbpd.lambda(end,:)); ...
    norm(cdbpd_t.eta(end,:) - cdbpd.eta(end,:)); ...
    abs(cdbpd_t.kkt(end) - cdbpd.kkt(end)); ...
    abs(cdbpd_t.ec(end) - cdbpd.ec(end)); ...
    abs(cdbpd_t.gap(end) - cdbpd.gap(end))];
solver_sensitivity = table(sens_qty, sens_diff, ...
    'VariableNames', {'Quantity','AbsDifference'});
fprintf('Solver sensitivity (baseline vs tightened): direct ||dx||=%.3e, C-DBPD ||dx||=%.3e\n', ...
    sens_diff(1), sens_diff(5));
disp(solver_sensitivity);
writetable(solver_sensitivity, fullfile(out_dir, 'nonglobal_solver_sensitivity.csv'));

assert(direct.kkt(end) < 1e-6 && norm(direct.x(end,:).' - d.x_local) < 1e-6 && direct.gap(end) > 100, ...
    'Direct DBPD did not converge accurately to the nonglobal KKT branch.');
assert(cdbpd.kkt(end) < 1e-6 && norm(cdbpd.x(end,:).' - d.x_global) < 1e-6 && ...
    cdbpd.ec(end) < 1e-6 && cdbpd.gap(end) < 1e-6, ...
    'Proposed C-DBPD did not converge accurately to the global vGNE.');
assert(direct_t.kkt(end) < 1e-6 && norm(direct_t.x(end,:).' - d.x_local) < 1e-6 && direct_t.gap(end) > 100, ...
    'Tightened Direct DBPD check failed.');
assert(cdbpd_t.kkt(end) < 1e-6 && norm(cdbpd_t.x(end,:).' - d.x_global) < 1e-6 && ...
    cdbpd_t.ec(end) < 1e-6 && cdbpd_t.gap(end) < 1e-6, ...
    'Tightened C-DBPD check failed.');

summary = table( ...
    ["Direct DBPD (5)"; "Proposed C-DBPD (11)"], ...
    [norm(direct.x(end,:).' - d.x_local); norm(cdbpd.x(end,:).' - d.x_global)], ...
    [direct.kkt(end); cdbpd.kkt(end)], ...
    [direct.orig_g(end); cdbpd.orig_g(end)], ...
    [NaN; cdbpd.ec(end)], ...
    [direct.gap(end); cdbpd.gap(end)], ...
    [max(abs(sum(direct.omega, 2))); max(abs(sum(cdbpd.omega, 2)))], ...
    'VariableNames', {'Algorithm','DistToBranch','TerminalKKT', ...
    'TerminalConstraint','TerminalConsistency','TerminalGlobalBRGap','MaxZeroSumOmega'});
disp(summary);
writetable(summary, fullfile(out_dir, 'nonglobal_branch_summary.csv'));
write_report(d, init, direct, cdbpd, solver_sensitivity, out_dir);
make_figures(direct, cdbpd, d, out_dir);

save(fullfile(out_dir, 'nonglobal_run_data.mat'), ...
    'd', 'init', 'direct', 'cdbpd', 'summary', 'solver_sensitivity', '-v7.3');
fprintf('Completed. Results: %s\n', out_dir);

%% Local functions
function d = build_model()
    d.N = 4; d.n = 1; d.d = 1;
    % Undirected 4-cycle: edges (1,2), (2,3), (3,4), (4,1)
    d.A = [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0];
    d.L = diag(sum(d.A, 2)) - d.A;
    ev = sort(real(eig(d.L))); d.lambda2 = ev(2); d.normL = ev(end);

    d.r_out = [1.32; 1.42; 1.50; 1.36];
    d.r_in  = [0.74; 0.82; 0.88; 0.78];
    d.w     = [8.0; 9.0; 8.5; 9.5];
    d.beta  = 1.0;
    d.lam_global = 3.0;
    d.lam_local  = 9.0;
    d.lam_bar    = 12.0;
    d.xlo = -2.5; d.xhi = 2.5;
    d.etalo = 1.0; d.etahi = 11.0;
    d.gamma = 25.0;
    d.mu_theta = 1.0; d.mu_varphi = 1.0; d.mu_phi = 1.0;

    d.a  = 0.5*(d.r_out.^2 + d.r_in.^2);
    d.dr = 0.5*(d.r_out.^2 - d.r_in.^2);
    d.B  = 0.5*d.w.*(d.dr.^2);

    rhs_q = 2*d.w.*d.dr.*(d.lam_local*d.r_in - d.lam_global*d.r_out) ...
        - d.beta*(d.L*(d.r_out + d.r_in));
    d.q = rhs_q ./ (d.r_out + d.r_in);
    assert(all(d.q > 10));
    d.r = d.r_out + (d.beta*(d.L*d.r_out) + ...
        d.lam_global*2*d.w.*d.dr.*d.r_out) ./ d.q;

    d.x_global = d.r_out;
    d.x_local  = -d.r_in;
    d.eta_star = d.w .* d.dr;
    assert(all(d.eta_star > d.etalo & d.eta_star < d.etahi));

    d.Romega = [eye(d.N-1); -ones(1, d.N-1)];
    d.Romega_pinv = (d.Romega.'*d.Romega)\d.Romega.';

    % Uniform canonical Slater at xhat = 0 for all eta_i in [etalo, etahi]
    d.slater_margin = -sum(d.etalo*(-d.a) - d.etalo^2./(2*d.w) - d.B);
    assert(d.slater_margin > 10);

    d.vi_at_local = F_map(d.x_local, d).' * (d.x_global - d.x_local);
    d.gap_at_local = global_br_gap(d.x_local, d);
    d.gap_at_global = global_br_gap(d.x_global, d);
    d.kkt_at_local = orig_kkt_residual(d.x_local, d.lam_local*ones(d.N,1), d);
    d.kkt_at_global = orig_kkt_residual(d.x_global, d.lam_global*ones(d.N,1), d);
end

function init = build_initialization(d)
    N = d.N;
    % Identical (x0, X0, lambda0, omega0) near the nonglobal KKT branch
    x0 = d.x_local + [0.05; -0.04; 0.04; -0.05];
    X0 = repmat(x0, 1, N) + 0.04*[ ...
        0.0,  0.5, -0.4,  0.3; ...
       -0.5,  0.0,  0.4, -0.3; ...
        0.3, -0.4,  0.0,  0.5; ...
       -0.3,  0.4, -0.5,  0.0];
    for i = 1:N, X0(i, i) = x0(i); end
    lambda0 = [8.80; 8.55; 8.75; 8.60];
    omega0 = zeros(N, 1);
    eta0 = [4.0; 5.0; 5.5; 4.8];

    Q0 = X0;
    for i = 1:N
        Q0(i, i) = fermi_inverse(x0(i), d.xlo, d.xhi, d.mu_theta);
    end
    nu0 = fermi_inverse(lambda0, 0, d.lam_bar, d.mu_varphi);
    zeta0 = omega0(1:N-1);
    alpha0 = fermi_inverse(eta0, d.etalo, d.etahi, d.mu_phi);

    init.x0 = x0; init.X0 = X0; init.lambda0 = lambda0;
    init.omega0 = omega0; init.eta0 = eta0;
    init.y0_direct = [Q0(:); nu0; zeta0];
    init.y0_cdbpd  = [Q0(:); nu0; zeta0; alpha0];
end

function dy = rhs_direct_dbpd(y, d)
    N = d.N; nq = N*N; nw = N - 1;
    Q = reshape(y(1:nq), N, N);
    nu = y(nq+(1:N));
    zeta = y(nq+N+(1:nw));
    om = d.Romega * zeta;
    X = Q; x = zeros(N, 1);
    for i = 1:N
        x(i) = fermi_map(Q(i, i), d.xlo, d.xhi, d.mu_theta);
        X(i, i) = x(i);
    end
    lam = fermi_map(nu, 0, d.lam_bar, d.mu_varphi);
    dQ = zeros(N, N); dnu = zeros(N, 1);
    for i = 1:N
        nbrs = find(d.A(i, :));
        Fi = d.q(i)*(x(i) - d.r(i)) + d.beta*sum(x(i) - X(nbrs, i));
        grad_gi = 2*d.w(i)*(x(i)^2 - d.a(i))*x(i);
        dQ(i, i) = -Fi - grad_gi*lam(i) - d.gamma*sum(x(i) - X(i, nbrs));
        for pl = 1:N
            if pl == i, continue; end
            dQ(pl, i) = -d.gamma*sum(X(pl, i) - X(pl, nbrs));
        end
        gi = 0.5*d.w(i)*(x(i)^2 - d.a(i))^2 - d.B(i);
        dnu(i) = gi - d.L(i,:)*om;
    end
    dzeta = d.Romega_pinv * (d.L * lam);
    dy = [dQ(:); dnu; dzeta];
end

function dy = rhs_cdbpd(y, d)
    N = d.N; nq = N*N; nw = N - 1;
    Q = reshape(y(1:nq), N, N);
    nu = y(nq+(1:N));
    zeta = y(nq+N+(1:nw));
    al = y(nq+N+nw+(1:N));
    om = d.Romega * zeta;
    X = Q; x = zeros(N, 1);
    for i = 1:N
        x(i) = fermi_map(Q(i, i), d.xlo, d.xhi, d.mu_theta);
        X(i, i) = x(i);
    end
    lam = fermi_map(nu, 0, d.lam_bar, d.mu_varphi);
    eta = fermi_map(al, d.etalo, d.etahi, d.mu_phi);
    dQ = zeros(N, N); dnu = zeros(N, 1); dal = zeros(N, 1);
    for i = 1:N
        nbrs = find(d.A(i, :));
        Fi = d.q(i)*(x(i) - d.r(i)) + d.beta*sum(x(i) - X(nbrs, i));
        dQ(i, i) = -Fi - 2*x(i)*eta(i)*lam(i) - d.gamma*sum(x(i) - X(i, nbrs));
        for pl = 1:N
            if pl == i, continue; end
            dQ(pl, i) = -d.gamma*sum(X(pl, i) - X(pl, nbrs));
        end
        rhoi = x(i)^2 - d.a(i);
        phistar = eta(i)^2/(2*d.w(i)) + d.B(i);
        dnu(i) = eta(i)*rhoi - phistar - d.L(i,:)*om;
        dal(i) = lam(i)*(rhoi - eta(i)/d.w(i));
    end
    dzeta = d.Romega_pinv * (d.L * lam);
    dy = [dQ(:); dnu; dzeta; dal];
end

function tr = decode_direct(t, y, d)
    nt = numel(t); N = d.N; nq = N*N; nw = N - 1;
    tr.t = t; tr.x = zeros(nt, N); tr.lambda = zeros(nt, N);
    tr.omega = zeros(nt, N); tr.kkt = zeros(nt, 1);
    tr.gap = zeros(nt, 1); tr.orig_g = zeros(nt, 1);
    for k = 1:nt
        Q = reshape(y(k, 1:nq), N, N);
        x = zeros(N, 1);
        for i = 1:N, x(i) = fermi_map(Q(i,i), d.xlo, d.xhi, d.mu_theta); end
        lam = fermi_map(y(k, nq+(1:N)).', 0, d.lam_bar, d.mu_varphi);
        om = d.Romega * y(k, nq+N+(1:nw)).';
        tr.x(k,:) = x.'; tr.lambda(k,:) = lam.'; tr.omega(k,:) = om.';
        tr.kkt(k) = orig_kkt_residual(x, lam, d);
        tr.gap(k) = global_br_gap(x, d);
        tr.orig_g(k) = sum(orig_g(x, d));
    end
end

function tr = decode_cdbpd(t, y, d)
    nt = numel(t); N = d.N; nq = N*N; nw = N - 1;
    tr.t = t; tr.x = zeros(nt, N); tr.lambda = zeros(nt, N);
    tr.omega = zeros(nt, N); tr.eta = zeros(nt, N);
    tr.kkt = zeros(nt, 1); tr.gap = zeros(nt, 1);
    tr.orig_g = zeros(nt, 1); tr.ec = zeros(nt, 1);
    for k = 1:nt
        Q = reshape(y(k, 1:nq), N, N);
        x = zeros(N, 1);
        for i = 1:N, x(i) = fermi_map(Q(i,i), d.xlo, d.xhi, d.mu_theta); end
        lam = fermi_map(y(k, nq+(1:N)).', 0, d.lam_bar, d.mu_varphi);
        om = d.Romega * y(k, nq+N+(1:nw)).';
        eta = fermi_map(y(k, nq+N+nw+(1:N)).', d.etalo, d.etahi, d.mu_phi);
        tr.x(k,:) = x.'; tr.lambda(k,:) = lam.'; tr.omega(k,:) = om.'; tr.eta(k,:) = eta.';
        tr.kkt(k) = orig_kkt_residual(x, lam, d);
        tr.gap(k) = global_br_gap(x, d);
        tr.orig_g(k) = sum(orig_g(x, d));
        tr.ec(k) = norm(eta - d.w.*(x.^2 - d.a));
    end
end

function F = F_map(x, d)
    F = d.q.*(x - d.r) + d.beta*(d.L*x);
end

function g = orig_g(x, d)
    g = 0.5*d.w.*(x.^2 - d.a).^2 - d.B;
end

function r = orig_kkt_residual(x, lam, d)
    gradg = 2*d.w.*(x.^2 - d.a).*x;
    G = sum(orig_g(x, d));
    lam_bar = mean(lam);
    r = norm(F_map(x, d) + lam.*gradg) + norm(lam - lam_bar) + abs(G);
end

function gap = global_br_gap(x, d)
    gaps = zeros(d.N, 1);
    g_all = orig_g(x, d);
    for i = 1:d.N
        nbrs = find(d.A(i, :));
        cur = 0.5*d.q(i)*(x(i) - d.r(i))^2 + 0.5*d.beta*sum((x(i) - x(nbrs)).^2);
        rhs = d.B(i) - (sum(g_all) - g_all(i));
        if rhs < 0
            gaps(i) = NaN; continue;
        end
        D = sqrt(max(0, 2*rhs/d.w(i)));
        rlo = sqrt(max(0, d.a(i) - D));
        rhi = min(d.xhi, sqrt(max(0, d.a(i) + D)));
        yfree = (d.q(i)*d.r(i) + d.beta*sum(x(nbrs))) / (d.q(i) + d.beta*numel(nbrs));
        cands = [-rhi; -rlo; rlo; rhi; ...
            min(max(yfree, -rhi), -rlo); min(max(yfree, rlo), rhi)];
        vals = zeros(size(cands));
        for m = 1:numel(cands)
            vals(m) = 0.5*d.q(i)*(cands(m) - d.r(i))^2 + ...
                0.5*d.beta*sum((cands(m) - x(nbrs)).^2);
        end
        gaps(i) = max(cur - min(vals), 0);
    end
    gap = max(gaps, [], 'omitnan');
end

function write_report(d, init, direct, cdbpd, sens, out_dir)
    fid = fopen(fullfile(out_dir, 'nonglobal_branch_report.txt'), 'w');
    assert(fid ~= -1);
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, 'EXAMPLE 3: DISTRIBUTED GLOBAL RECOVERY VS NONGLOBAL KKT BRANCH\n');
    fprintf(fid, 'MATLAB version: %s\n', version);
    fprintf(fid, 'Graph: 4-cycle, edges (1,2),(2,3),(3,4),(4,1)\n');
    fprintf(fid, 'Solver: ode15s, RelTol=2e-10, AbsTol=2e-12, MaxStep=0.02, T=80\n');
    fprintf(fid, 'q = %.17g %.17g %.17g %.17g\n', d.q);
    fprintf(fid, 'r = %.17g %.17g %.17g %.17g\n', d.r);
    fprintf(fid, 'Identical x0 = %.17g %.17g %.17g %.17g\n', init.x0);
    fprintf(fid, 'Identical lambda0 = %.17g %.17g %.17g %.17g\n', init.lambda0);
    fprintf(fid, 'Identical omega0 = %.17g %.17g %.17g %.17g\n', init.omega0);
    fprintf(fid, 'C-DBPD eta0 = %.17g %.17g %.17g %.17g\n', init.eta0);
    fprintf(fid, 'Direct DBPD terminal x = %.17g %.17g %.17g %.17g\n', direct.x(end,:));
    fprintf(fid, 'Direct DBPD terminal lambda = %.17g %.17g %.17g %.17g\n', direct.lambda(end,:));
    fprintf(fid, 'Direct DBPD terminal KKT = %.17g, G_BR = %.17g\n', direct.kkt(end), direct.gap(end));
    fprintf(fid, 'C-DBPD terminal x = %.17g %.17g %.17g %.17g\n', cdbpd.x(end,:));
    fprintf(fid, 'C-DBPD terminal lambda = %.17g %.17g %.17g %.17g\n', cdbpd.lambda(end,:));
    fprintf(fid, 'C-DBPD terminal eta = %.17g %.17g %.17g %.17g\n', cdbpd.eta(end,:));
    fprintf(fid, 'C-DBPD terminal KKT = %.17g, e_c = %.17g, G_BR = %.17g\n', ...
        cdbpd.kkt(end), cdbpd.ec(end), cdbpd.gap(end));
    fprintf(fid, 'Tightened-solver max |dx|: direct = %.17g, C-DBPD = %.17g\n', ...
        sens.AbsDifference(1), sens.AbsDifference(5));
end

function make_figures(direct, cdbpd, d, out_dir)
    % Exact first 4 node colors matching the convex figures
    cols = [ ...
        0.0000, 0.0000, 1.0000; ...
        1.0000, 0.0000, 0.0000; ...
        0.0000, 0.6275, 0.0000; ...
        1.0000, 0.5490, 0.0000];

    f = figure('Color', 'w', 'Units', 'pixels', 'Position', [1 1 2048 1030], ...
        'Visible', 'off');
    set(f, 'PaperUnits', 'points', 'PaperPosition', [0 0 1485 750], ...
        'PaperSize', [1485 750]);
    tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Left: strategy trajectories of all 4 agents
    ax = nexttile(tl); hold(ax, 'on');
    hc = gobjects(d.N, 1); hd = gobjects(d.N, 1);
    for i = 1:d.N
        hc(i) = plot(ax, cdbpd.t, cdbpd.x(:, i), 'LineStyle', '-', ...
            'Color', cols(i,:), 'LineWidth', 5.0, ...
            'DisplayName', sprintf('C-DBPD: $x_%d$', i));
        hd(i) = plot(ax, direct.t, direct.x(:, i), 'LineStyle', '-.', ...
            'Color', cols(i,:), 'LineWidth', 5.0, ...
            'DisplayName', sprintf('DBPD: $x_%d$', i));
    end
    yline(ax, 0, 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2.5, ...
        'HandleVisibility', 'off');
    xlim(ax, [0, 10]); ylim(ax, [-1.05, 1.65]);
    xticks(ax, [0, 5, 10]); yticks(ax, [-1, -0.5, 0, 0.5, 1, 1.5]);
    xlabel(ax, 'Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', 55);
    ylabel(ax, '$x_i$', 'Interpreter', 'latex', 'FontSize', 55);
    lgd = legend(ax, [hc; hd], 'Interpreter', 'latex', ...
        'FontSize', 33, 'NumColumns', 2, 'Box', 'off');
    lgd.ItemTokenSize = [24, 12];
    lgd.Position = [0.148, 0.495, 0.320, 0.300];

    % Right: unilateral global best-response gap and KKT residual
    ax = nexttile(tl); hold(ax, 'on');
    h1 = plot(ax, direct.t, max(direct.gap, 1e-12), 'LineStyle', '--', ...
        'LineWidth', 5.0, 'Color', [0.0000 0.0000 1.0000], ...
        'DisplayName', 'DBPD: $\mathcal{G}_{\mathrm{BR}}$');
    h2 = plot(ax, cdbpd.t, max(cdbpd.gap, 1e-12), 'LineStyle', '-', ...
        'LineWidth', 5.0, 'Color', [1.0000 0.0000 0.0000], ...
        'DisplayName', 'C-DBPD: $\mathcal{G}_{\mathrm{BR}}$');
    h3 = plot(ax, direct.t, max(direct.kkt, 1e-12), 'LineStyle', ':', ...
        'LineWidth', 5.0, 'Color', [0.0000 0.7216 0.8510], ...
        'DisplayName', 'DBPD: KKT');
    h4 = plot(ax, cdbpd.t, max(cdbpd.kkt, 1e-12), 'LineStyle', '-.', ...
        'LineWidth', 4.5, 'Color', [0.0000 0.6275 0.0000], ...
        'DisplayName', 'C-DBPD: KKT');
    set(ax, 'YScale', 'log');
    xlim(ax, [0, 80]); ylim(ax, [1e-10, 1e3]);
    xticks(ax, [0, 20, 40, 60, 80]);
    yticks(ax, [1e-10, 1e-6, 1e-2, 1e2]);
    xlabel(ax, 'Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', 55);
    ylabel(ax, '$\mathcal{G}_{\mathrm{BR}}(x(t))$ and KKT', ...
        'Interpreter', 'latex', 'FontSize', 55);
    lgd2 = legend(ax, [h1 h2 h3 h4], 'Location', 'southwest', ...
        'Interpreter', 'latex', 'FontSize', 42, 'Box', 'off');
    lgd2.ItemTokenSize = [28, 16];

    ax_all = findall(f, 'Type', 'axes');
    set(ax_all, 'FontName', 'Times New Roman', 'FontSize', 50, 'LineWidth', 0.90, ...
        'TickLabelInterpreter', 'latex', 'Box', 'on', ...
        'XGrid', 'off', 'YGrid', 'off', 'YMinorGrid', 'off');

    set(f, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    savefig(f, fullfile(out_dir, 'nonconvex_globality.fig'));
    exportgraphics(f, fullfile(out_dir, 'nonconvex_globality.png'), 'Resolution', 200);
    try
        exportgraphics(f, fullfile(out_dir, 'nonconvex_globality.eps'), 'ContentType', 'vector');
    catch
        print(f, fullfile(out_dir, 'nonconvex_globality.eps'), '-depsc', '-painters');
    end
    close(f);
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

function v = expand(v, n)
    if isscalar(v), v = repmat(v, n, 1); else, v = v(:); end
end
