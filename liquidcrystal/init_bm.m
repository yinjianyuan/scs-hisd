%% precompute initial states and eigenvectors once for fairness
init_dir = './init_benchmark';
if ~exist(init_dir, 'dir')
    mkdir(init_dir);
end

for ic = 1:numel(cases)
    tag = cases(ic).tag;
    init_file = fullfile(init_dir, ['init_' tag '.mat']);

    if exist(init_file, 'file')
        fprintf('Loading fixed initialization for %s from %s...\n', tag, init_file);
        Sinit = load(init_file, 'u', 'Vfull', 'alpha', 'x0');

        cases(ic).u = Sinit.u;
        cases(ic).Vfull = Sinit.Vfull;
        cases(ic).alpha = Sinit.alpha;
        cases(ic).x0 = Sinit.x0;
    else
        error("No such file.")
    end
end