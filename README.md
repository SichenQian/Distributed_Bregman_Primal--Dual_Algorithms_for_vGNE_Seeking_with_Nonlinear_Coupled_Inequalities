# Distributed Bregman Primal–Dual Algorithms for vGNE Seeking with Nonlinear Coupled Inequalities

Python and MATLAB implementations and reproducibility package for the paper:

**Distributed Bregman Primal–Dual Algorithms for vGNE Seeking with Nonlinear Coupled Inequalities**  
Sichen Qian, Hongzhe Liu, Wenwu Yu, and Wei Xing Zheng.

**Paper status:** Prepared for submission to *IEEE Transactions on Automatic Control*.

## Numerical examples

| Script | Example | Manuscript items |
| --- | --- | --- |
| `main.py` | Two-player convex and structured composite nonconvex feasible sets, potential contours, and equilibrium points | Fig. 1 |
| `example_1_convex.m` | Ten-agent edge-server regulation game with a convex nonlinear coupled inequality | Figs. 2–5, Table II |
| `example2_nonconvex.m` | Ten-agent edge-server regulation game with a structured composite nonconvex coupled inequality | Figs. 6–8, Table III |
| `example3_nonglobal_branch.m` | Four-agent coupled game over a four-cycle contrasting direct DBPD and C-DBPD from identical initialization near a nonglobal KKT branch | Fig. 9 |

All three MATLAB scripts are self-contained. Model parameters, communication graphs, equilibrium targets, and deterministic initializations are defined inside each script. Required theorem inequalities and post-simulation residual checks are enforced via assertions.

## Requirements

### Python

- Python **3.12 or later**; validated with **3.13.15**.
- NumPy and Matplotlib, with versions specified in `requirements.txt`.

The default Python run uses Matplotlib's built-in mathematical text rendering and bundled STIX fonts, so no external LaTeX installation is required.

### MATLAB

- MATLAB; validated with **R2024a**.
- No additional MATLAB toolboxes, Python packages, or external input files are required.

## Usage

### Python: Fig. 1

From the repository folder, install the dependencies and run:

```sh
python -m pip install -r requirements.txt
python main.py
```

The script saves the two panels as `convex_case.eps`, `convex_case.png`, `nonconvex_case.eps`, and `nonconvex_case.png` in `results/figure1/`.

An alternative output folder can be selected with `--output-dir`:

```sh
python main.py --output-dir path/to/results
```

For the manuscript's LaTeX font configuration, use:

```sh
python main.py --usetex
```

This optional mode requires `latex`, `dvips`, and `dvipng` on `PATH`, together with the LaTeX packages `amsmath`, `amssymb`, `newtxtext`, and `newtxmath`. Ghostscript is used for EPS distillation when available. See the [Matplotlib LaTeX rendering documentation](https://matplotlib.org/stable/users/explain/text/usetex.html).

### MATLAB: Figs. 2–9

Open any of the MATLAB scripts and click **Run**, or execute the following commands from the repository folder:

```matlab
example_1_convex
example2_nonconvex
example3_nonglobal_branch
```

All four programs save their default outputs relative to their own location, so no output paths need to be edited:

```text
results/
  figure1/
  example1_convex/
  example2_nonconvex/
  example3_nonglobal/
```

Running a program again overwrites only its corresponding output folder.

## Outputs

Each MATLAB script exports figures in `.fig`, `.png`, and `.eps` formats, CSV audit tables, and a MAT file containing model parameters, initialization, trajectories, and numerical metrics:

- `results/example1_convex/convex_run_data.mat`
- `results/example2_nonconvex/nonconvex_run_data.mat`
- `results/example3_nonglobal/nonglobal_run_data.mat`

Both `example2_nonconvex.m` and `example3_nonglobal_branch.m` also execute a second integration under tightened solver tolerances and a halved maximum step size, saving solver-sensitivity comparison tables alongside the baseline outputs. Every script terminates with an assertion error if any required check fails.

## Figure correspondence

The Python panels are exported as `.eps` and `.png`; the MATLAB figures are exported as `.fig`, `.png`, and `.eps`:

| Figure | Output folder | File stem | Content |
| --- | --- | --- | --- |
| Fig. 1a | `results/figure1/` | `convex_case` | Convex feasible set, potential contours, and vGNE |
| Fig. 1b | `results/figure1/` | `nonconvex_case` | Nonconvex feasible set, potential contours, global vGNE, and nonglobal KKT point |
| Fig. 2 | `results/example1_convex/` | `shared_topology` | Shared ten-node communication graph |
| Fig. 3 | `results/example1_convex/` | `convex_tra` | Primal and decision-estimate squared errors |
| Fig. 4 | `results/example1_convex/` | `convex_aux` | Multiplier consensus and auxiliary-state trajectories |
| Fig. 5 | `results/example1_convex/` | `convex_exp` | Normalized logarithmic aggregate error and shared-constraint residual |
| Fig. 6 | `results/example2_nonconvex/` | `nonconvex_tra` | Primal and decision-estimate squared errors for the ten-agent nonconvex example |
| Fig. 7 | `results/example2_nonconvex/` | `nonconvex_aux` | Scaled multiplier deviations and scaled auxiliary states |
| Fig. 8 | `results/example2_nonconvex/` | `nonconvex_exp` | Normalized logarithmic aggregate, canonical-state, and consistency errors; original and canonical constraint residuals |
| Fig. 9 | `results/example3_nonglobal/` | `nonconvex_globality` | Four-agent strategy trajectories, unilateral global best-response gaps, and KKT residuals |

## Numerical verification

### Certificate verification

The numerical certificates for Tables II and III are generated by `example_1_convex.m` and `example2_nonconvex.m`, respectively:

| Manuscript item | Verification output |
| --- | --- |
| Table II | `results/example1_convex/convex_theorem_audit.csv` |
| Table III | `results/example2_nonconvex/nonconvex_theorem_audit.csv` |
| Table III full precision | `results/example2_nonconvex/table_III_certificate.txt` |

The Table III text certificate records 17 significant digits of the double-precision values used to verify every strict inequality in Theorem 2 and Appendix G, including the core small-gain condition, the auxiliary parameter bounds, the selected gain `gamma = 2200`, the rank margin `iota = 2025`, and all five Lyapunov positivity margins. Infinite one-sided bounds appear as `Inf` or `-Inf`.

The script evaluates the certificate using computed model and sublevel constants stored in `results/example2_nonconvex/nonconvex_model_constants.csv`, whereas the manuscript reports conservative rounded upper bounds. The code-computed threshold for the rank condition in Eq. 31 is approximately `1.5287660333e3`; substituting the manuscript's displayed rounded upper bounds yields approximately `1.5287914450e3`. Both values lie below the displayed bound `1.529e3` and well below `iota = 2025`.

### Four-agent global-recovery comparison

Fig. 9 is generated by `example3_nonglobal_branch.m`, which simulates the direct DBPD algorithm on the original nonconvex constraint and the proposed C-DBPD algorithm over an undirected four-cycle from identical primal, opponent-estimate, multiplier, and auxiliary-state initializations near a locally stable nonglobal KKT branch. This comparison is separate from the ten-agent regional certificate in `example2_nonconvex.m` and illustrates the global-recovery mechanism of Lemma 6.

Full-precision model parameters, initial states, terminal states, KKT residuals, consistency residuals, unilateral global best-response gaps, and tightened-solver differences are recorded in:

- `results/example3_nonglobal/nonglobal_branch_report.txt`
- `results/example3_nonglobal/nonglobal_branch_summary.csv`
- `results/example3_nonglobal/nonglobal_solver_sensitivity.csv`

### Verification records

From the repository folder, record a full MATLAB run with a fresh diary filename so MATLAB does not append to an existing log:

```matlab
diary('reproduction_log_local.txt')
disp('Reproducing Tables II-III and Figs. 2-9')
disp(version)
example_1_convex
example2_nonconvex
example3_nonglobal_branch
diary off
```

Run `python main.py` separately for Fig. 1 after installing `requirements.txt`. The release includes `reproduction_log.txt`, `python_reproduction_log.txt`, all theorem and post-simulation audit tables, and the solver-sensitivity tables. All three MATLAB scripts must reach `Completed. Results:`; any failed theorem, baseline, tightened-solver, or branch-recovery assertion stops execution.

The fixed release is [`TAC-v9-reproducibility`](https://github.com/SichenQian/Distributed_Bregman_Primal--Dual_Algorithms_for_vGNE_Seeking_with_Nonlinear_Coupled_Inequalities/releases/tag/TAC-v9-reproducibility). The verification logs supplied with this release correspond to the tagged source-code version. `verification_manifest.json` records the runtime versions, validation counts, high-precision bound comparisons, and Git blob hashes of the source files, documentation, and verification records. To check out the fixed source version:

```sh
git checkout TAC-v9-reproducibility
git rev-parse HEAD
```

Generated figures and large MAT files are excluded from Git and are regenerated by the commands above.

## Citation

If you use this code in your research, please cite the accompanying paper:

```bibtex
@unpublished{qian_bregman_vgne,
  author = {Qian, Sichen and Liu, Hongzhe and Yu, Wenwu and Zheng, Wei Xing},
  title  = {Distributed {Bregman} Primal--Dual Algorithms for {vGNE} Seeking with Nonlinear Coupled Inequalities},
  note   = {Submitted to {IEEE} Transactions on Automatic Control}
}
```
