H=@hv;

optionsini = hiosdoptions(1);
optionsini.dt = 0;       
optionsini.outputX = 0; 
optionsini.outputp = 100;  
optionsini.outputd = 0;
optionsini.maxiter =300;  
optionsini.k=  15; 
rng(0);
optionsini.V=randn(d,optionsini.k);
optionsini.s=s1;

[VV, alpha] = hiosd_lobpcg_iniP_eh(@hv,x,optionsini); 
