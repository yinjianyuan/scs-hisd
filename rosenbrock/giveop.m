F = @enegrad;
H = @hv;

options = hiosdoptions(1);
options.dt = options.dt;
options.betat = 8;
options.betau = 0.5;
options.gamma = 0.75;
options.s = s1;
options.epsf = 1e-10; 
options.maxiter = 800;
options.outputp = 1000;
options.outputd = 1000;
options=rmfield(options,{'gammamax', 'gammamin', 'l', 'minl', 'tau'});