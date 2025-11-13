
clear all
close all
ymin = 2.0;
ymax = 2.025;
xmin = 39;
xmax = 40;
file = "MS_tauD_1_to_7_NtauD_20_sigmaD_m0p4_to_2_NsigmaD_80.mat";

[sigma,tau,psi_L] = Extract_values(xmin,xmax,ymin,ymax,file)
clc

sigma
tau
psi_L
