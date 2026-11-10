clear;clc
initialize 
initeigs;
giveop;
options.V = VV(:,1:5);
[x1,output1] = hisd_lobpcg_nobb(F, H, x, options);
strictly_euler1 = output1.data;

eta0=2.6;
[x1,output2] = scshisd_lobpcg_nobb(F, H, x, options, eta0);
strictly_euler2 = output2.data;

draw;