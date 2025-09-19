clear all; close all; clc; 
% script_test_pseudo_resonant_bug.m

 load('Case_81_transfer_bugData.mat')
% load('Case_18_transfer_bugData.mat')

        N=2
        M=3

        vinf=TitanNODE_Arrival(1);
        idmoon=5;
        idcentral=6;

        [vinf1, alpha1, crank1, vinf2, alpha2, crank2, tofsc] = wrap_pseudoResTransf(TypeFlag, N, M, vinf, idmoon, idcentral);

        [~, rr, vv] = vinfAlphaCrank2car([vinf1 alpha1 crank1], SVDeparture, pars.Planet.mu);
        [kep_sc_out] = car2kep([rr vv], pars.Planet.mu);
        Period_spacecraft_out = 2*pi*sqrt(kep_sc_out(1)^3/pars.Planet.mu);
        
        tofsc/Period_titan
        tofsc/Period_spacecraft_out


        % now the bug seems repaired.