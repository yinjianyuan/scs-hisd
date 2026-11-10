
F = @enegrad;
H = @hv;

options = hiosdoptions(1);
options.dt = 0.02;
options.betat = nan;
options.betau = nan;
options.s = 0;
options.epsf = 1e-10; 
options.maxiter = 100000000;
options.outputp = 500;
options.outputd = 1000;
options=rmfield(options,{'gammamax', 'gammamin', 'l', 'minl', 'tau'});