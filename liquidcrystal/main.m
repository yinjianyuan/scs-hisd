clear; clc;
global I N s omega x y h A B C E0 A2C A4C B2C h2 xi
h = 1/50;  x = -1 : h : 1;   y = -1 : h : 1;   I = length(x) - 2;   N = I ^ 2;
%parameters
B = 0.64 * 10^4;  C = 0.35 * 10^4;   A = -B^2/3/C;   
s = ( B + sqrt(B^2 + 24 * abs(A) * C) )/4/C; xi = -B^2/2/C/C;
E0=s^2*A/(6*C)-s^3*B/(27*C)+s^4/18;
A2C=A/2/C; A4C=A/4/C; B2C=B/2/C; h2=h^2;
G=[2*ones(N,1);2*ones(N,1)];
inpf=@(x,y)x'*(G.*y)*h2;

omega =     50   ;

F=@(uu)ngradnlc_2d(uu);   

options=hiosdoptions(1);options.outputX=100; options.dt=0.0008;
options.gammamax=0.1;options.gammamin=0.001;options.betat=5;options.betau=0.2;
options.innpfunc=inpf;   options.epsf=1e-6;


%% 3-s
% % NAME='T1';
nRep = 1;                 % repeated runs
maxiter_bench = 8000000;   
force_full_iter = false;  

methods = {'SCSHiSD','AHiSD', 'HiSD'};

% Rotate method order to reduce order/JIT/cache bias
method_orders = {
    {'SCSHiSD','AHiSD', 'HiSD'}, ...
};

%% timing-related options
if exist('options','var')
    options_template = options;
else
    options_template = struct();
end

options_template.outputX = 0;          % avoid storing trajectories during timing
options_template.outputd = 0;          % avoid drawing
options_template.outputp = 10^12;      % avoid printing
options_template.maxiter = maxiter_bench;

if force_full_iter
    options_template.epsf = 0;          % norm is nonnegative, so this prevents early stop
end
disp(options_template.epsf)

% BB bounds: beta in [beta0/5, 5 beta0]
options_template.betat = 5;
options_template.betau = 1/5;

%% define two saddle searches
cases = struct();

% M+ -> N+
cases(1).tag = 'M_to_N';
cases(1).state_name = 'M+e';
cases(1).perturb_idx = 4;
cases(1).scs_k = 6;       % SCSHiSD+l, here k+l=6
cases(1).ind = 4;
cases(1).hisd_k = 4;      % HiSD/A-HiSD index
cases(1).dt_hisd = 2.4e-4;
cases(1).dt_scs = 2.4e-4;
cases(1).dt_ahisd = 8.0e-4;

% N+ -> M+
cases(2).tag = 'N_to_M';
cases(2).state_name = 'N+1';
cases(2).perturb_idx = 4;
cases(2).scs_k = 6;       % SCSHiSD+l, here k+l=6
cases(2).ind = 3;
cases(2).hisd_k = 3;      % according to the previous table: N+ -> M+ uses k=3
cases(2).dt_hisd = 2.4e-4;
cases(2).dt_scs = 2.4e-4;
cases(2).dt_ahisd = 8.0e-4;

cases(3).tag = 'T+_to_T';
cases(3).state_name = 'T+1';
cases(3).perturb_idx = 4;
cases(3).scs_k = 6;       % SCSHiSD+l, here k+l=6
cases(3).ind = 4;
cases(3).hisd_k = 4;      % according to the previous table: N+ -> M+ uses k=3
cases(3).dt_hisd = 2.4e-4;
cases(3).dt_scs = 2.4e-4;
cases(3).dt_ahisd = 8.0e-4;

cases(4).tag = 'T_to_T+';
cases(4).state_name = 'T1';
cases(4).perturb_idx = 4;
cases(4).scs_k = 6;       % SCSHiSD+l, here k+l=6
cases(4).ind = 3;
cases(4).hisd_k = 3;      % according to the previous table: N+ -> M+ uses k=3
cases(4).dt_hisd = 2.4e-4;
cases(4).dt_scs = 2.4e-4;
cases(4).dt_ahisd = 8.0e-4;

% M+ -> N+
cases(5).tag = 'C_to_M';
cases(5).state_name = 'C+';
cases(5).perturb_idx = 3;
cases(5).scs_k = 6;       % SCSHiSD+l, here k+l=6
cases(5).ind = 3;
cases(5).hisd_k = 3;      % HiSD/A-HiSD index
cases(5).dt_hisd = 2.4e-4;
cases(5).dt_scs = 2.4e-4;
cases(5).dt_ahisd = 8.0e-4;

% N+ -> M+
cases(6).tag = 'M_to_C';
cases(6).state_name = 'M+e';
cases(6).perturb_idx = 3;
cases(6).scs_k = 6;       % SCSHiSD+l, here k+l=6
cases(6).ind = 2;
cases(6).hisd_k = 2;      % according to the previous table: N+ -> M+ uses k=3
cases(6).dt_hisd = 2.4e-4;
cases(6).dt_scs = 2.4e-4;
cases(6).dt_ahisd = 8.0e-4;

%% initial benchmark
init_dir = './init_benchmark';
init_bm;

%% optional warm-up: avoids counting first-call JIT overhead
do_warmup = false; %%%%%
if do_warmup
    fprintf('\nWarm-up runs...\n');
    for ic = 1:numel(cases)
        for im = 1:numel(methods)
            opt = make_method_options(options_template, cases(ic), methods{im}, 20);
            xinit = cases(ic).x0;
            try
                run_one_method(F, xinit, opt, methods{im});
            catch ME
                fprintf('Warm-up failed for %s, %s: %s\n', cases(ic).tag, methods{im}, ME.message);
            end
        end
    end
end

%% timed runs
nTotal = nRep * numel(cases) * numel(methods);
Transition = cell(nTotal,1);
Method = cell(nTotal,1);
Rep = zeros(nTotal,1);
Iter = zeros(nTotal,1);
WallTime = zeros(nTotal,1);
CpuTime = zeros(nTotal,1);
WallPerIter = zeros(nTotal,1);
CpuPerIter = zeros(nTotal,1);

row = 0;

fprintf('\nTimed benchmark starts...\n');

for rep = 1:nRep
    order = method_orders{mod(rep-1,numel(method_orders))+1};

    fprintf('\n==================== Repeat %d / %d ====================\n', rep, nRep);

    for ic = 1:numel(cases)
        fprintf('\n----- Case %s -----\n', cases(ic).tag);

        for im = 1:numel(order)
            method = order{im};

            opt = make_method_options(options_template, cases(ic), method, maxiter_bench);
            xinit = cases(ic).x0;

            fprintf('Running %-8s | transition=%s | rep=%d ...\n', method, cases(ic).tag, rep);

            t0_wall = tic;
            t0_cpu = cputime;

            [~, out] = run_one_method(F, xinit, opt, method);

            wall_time = toc(t0_wall);
            cpu_time = cputime - t0_cpu;

            if isstruct(out) && isfield(out,'it')
                it = out.it;
            else
                it = maxiter_bench;
            end

            row = row + 1;
            Transition{row} = cases(ic).tag;
            Method{row} = method;
            Rep(row) = rep;
            Iter(row) = it;
            WallTime(row) = wall_time;
            CpuTime(row) = cpu_time;
            WallPerIter(row) = wall_time / max(it,1);
            CpuPerIter(row) = cpu_time / max(it,1);

            fprintf('%-8s | iter=%d | wall=%.4f s | cpu=%.4f s | cpu/iter=%.4e s\n', ...
                method, it, wall_time, cpu_time, CpuPerIter(row));

            benchmark_results = struct();
            benchmark_results.Transition = Transition(1:row);
            benchmark_results.Method = Method(1:row);
            benchmark_results.Rep = Rep(1:row);
            benchmark_results.Iter = Iter(1:row);
            benchmark_results.WallTime = WallTime(1:row);
            benchmark_results.CpuTime = CpuTime(1:row);
            benchmark_results.WallPerIter = WallPerIter(1:row);
            benchmark_results.CpuPerIter = CpuPerIter(1:row);
            save('benchmark_nlc_cpu_10000_partial.mat','benchmark_results');
        end
    end
end

%% collect table
T = table(Transition, Method, Rep, Iter, WallTime, CpuTime, WallPerIter, CpuPerIter);
writetable(T,'benchmark_nlc_cpu_10000_raw.csv');

%% summary table
Summary = make_summary_table(T, cases, methods);
disp(Summary);
writetable(Summary,'benchmark_nlc_cpu_10000_summary.csv');

save('benchmark_nlc_cpu_10000.mat','T','Summary','cases');

fprintf('\nSaved:\n');
fprintf('  benchmark_nlc_cpu_10000_raw.csv\n');
fprintf('  benchmark_nlc_cpu_10000_summary.csv\n');
fprintf('  benchmark_nlc_cpu_10000.mat\n');

%% helper functions
function opt = make_method_options(options_template, caseinfo, method, maxiter_now)
    opt = options_template;
    opt.maxiter = maxiter_now;
    opt.outputX = 0;
    opt.outputd = 0;
    opt.outputp = 10^12;
    opt.betat = 5;
    opt.betau = 1/5;
    opt.epsf = 1e-6;

    switch method
        case 'SCSHiSD'
            opt.k = caseinfo.scs_k;
            opt.ind = caseinfo.ind;
            opt.V = caseinfo.Vfull(:,1:caseinfo.scs_k);
            opt.dt = caseinfo.dt_scs;

        case 'HiSD'
            opt.k = caseinfo.hisd_k;
            opt.V = caseinfo.Vfull(:,1:caseinfo.hisd_k);
            opt.dt = caseinfo.dt_hisd;

        case 'AHiSD'
            opt.k = caseinfo.hisd_k;
            opt.V = caseinfo.Vfull(:,1:caseinfo.hisd_k);
            opt.dt = caseinfo.dt_ahisd;
            opt.gamma = 0.99;
            opt.gmma = 0.99;   % in case your local A-HiSD code uses gmma instead of gamma

        otherwise
            error('Unknown method: %s', method);
    end
end

function [u1, output] = run_one_method(F, xinit, opt, method)
    switch method
        case 'SCSHiSD'
            [u1, output] = hiosd_lobpcg_scs_posi(F, xinit, opt);

        case 'HiSD'
            [u1, output] = hiosd_lobpcg(F, xinit, opt);

        case 'AHiSD'
            [u1, output] = hisd_heavyball_bb(F, xinit, opt);

        otherwise
            error('Unknown method: %s', method);
    end
end

function Summary = make_summary_table(T, cases, methods)
    nrow = numel(cases) * numel(methods);

    Transition = cell(nrow,1);
    Method = cell(nrow,1);

    MeanCpu = zeros(nrow,1);
    MedianCpu = zeros(nrow,1);
    StdCpu = zeros(nrow,1);

    MeanWall = zeros(nrow,1);
    MedianWall = zeros(nrow,1);
    StdWall = zeros(nrow,1);

    MeanCpuPerIter = zeros(nrow,1);
    MedianCpuPerIter = zeros(nrow,1);
    MeanWallPerIter = zeros(nrow,1);
    MedianWallPerIter = zeros(nrow,1);

    MeanIter = zeros(nrow,1);

    row = 0;
    for ic = 1:numel(cases)
        for im = 1:numel(methods)
            row = row + 1;

            idx = strcmp(T.Transition, cases(ic).tag) & strcmp(T.Method, methods{im});

            Transition{row} = cases(ic).tag;
            Method{row} = methods{im};

            MeanCpu(row) = mean(T.CpuTime(idx));
            MedianCpu(row) = median(T.CpuTime(idx));
            StdCpu(row) = std(T.CpuTime(idx));

            MeanWall(row) = mean(T.WallTime(idx));
            MedianWall(row) = median(T.WallTime(idx));
            StdWall(row) = std(T.WallTime(idx));

            MeanCpuPerIter(row) = mean(T.CpuPerIter(idx));
            MedianCpuPerIter(row) = median(T.CpuPerIter(idx));

            MeanWallPerIter(row) = mean(T.WallPerIter(idx));
            MedianWallPerIter(row) = median(T.WallPerIter(idx));

            MeanIter(row) = mean(T.Iter(idx));
        end
    end

    Summary = table(Transition, Method, MeanIter, ...
        MeanCpu, MedianCpu, StdCpu, MeanCpuPerIter, MedianCpuPerIter, ...
        MeanWall, MedianWall, StdWall, MeanWallPerIter, MedianWallPerIter);
end
