load testbug.mat

N=1;
        M=1;

        vinf=TitanNODE_Arrival(1);
        idmoon=5;
        idcentral=6;

% --> solve the pseudo-resonant transfer
S                           = [ TypeFlag +1 1 1 0 ];
[~, ~, tofsc1, node1, alpha11] = wrap_VILT(S, vinf, vinf,  [], [], idmoon, idcentral);

S                           = [ TypeFlag +1 1 2 0 ];
[~, ~, tofsc2, node2, alpha12] = wrap_VILT(S, vinf, vinf,  [], [], idmoon, idcentral);


% It does not have any sense that the 1:1 resonance is longer than the 1:2
% resonance for this 81 transfer. 


        [vinf1, alpha1, crank1, vinf2, alpha2, crank2, tofsc] = wrap_pseudoResTransf(TypeFlag, N, M, vinf, idmoon, idcentral );

        [~, rr, vv] = vinfAlphaCrank2car([vinf1 alpha1 crank1], SVDeparture, pars.Planet.mu);
        [kep_sc_out] = car2kep([rr vv], pars.Planet.mu);
        Period_spacecraft_out = 2*pi*sqrt(kep_sc(1)^3/pars.Planet.mu);

     
        tofsc/Period_titan
        tofsc/Period_spacecraft_out