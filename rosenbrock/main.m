clear; clc; close all;
casenumber = 1;
s1s = [-662;-659;-657.5;-657.282];
initialize;
load("x.mat")
initeigs;
giveop;

options.V = VV(:,1:4);
options.rao = 1;

options.dt=0.001; eta0=180;
[x1,output1] = scshisd_nobb(F, H, x, options, eta0);
scs_nobb = output1.data;

options.dt=0.001;
[x1,output1] = hisd_nobb(F, H, x, options);
hisd_nobb = output1.data;

options.dt=0.001;
[x1,output1] = ahisd_nobb(F, H, x, options);
ahisd_nobb = output1.data;

options.dt=0.001;
[x1,output1] = scshisd_bb(F, H, x, options, eta0);
scs_bb = output1.data;

options.dt=0.001;
[x1,output1] = hisd_bb(F, H, x, options);
hisd_bb = output1.data;

options.dt=0.001;
options.betau = 1;
[x1,output1] = ahisd_bb(F, H, x, options);
ahisd_bb = output1.data;

draw;
