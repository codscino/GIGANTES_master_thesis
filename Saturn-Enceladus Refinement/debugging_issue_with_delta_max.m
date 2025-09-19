    %Determine maximum bending due to flyby
    rp_flyby  = pars.INPUTS.Flyby.min_h(1) + pars.Moon.EquRad(1);           %[km
    e_fly     = 1 + ((rp_flyby*vinfin^2)/pars.Moon.mu(1));    %[-]
    delta_max = 2*asin(1/e_fly);                                      %[rad]


% Converting Incoming Node into Planeto Centric Cartesian state vector in TCN frame
[vvinfin, rr_in, vv_in, ~] = vinfAlphaCrank_to_VinfCART(vinfin, alphain, kin, epoch0, pars.INPUTS.idMoon, pars.INPUTS.idCentral);
 

KepEnceladus=[pars.Moon.OrbRad(1) 0 0 0 0 0];
SVEnceladus = @(TrA)kep2car([KepEnceladus(1:5) TrA],pars.Planet.mu);
 SVArrival=SVEnceladus(0);
v_inf_minus_enceladus=vinfAlphaCrank2car(nodein, SVArrival, pars.Planet.mu);

v_inf_plus_enceladus=vinfAlphaCrank2car(nodeout, SVArrival, pars.Planet.mu);


delta = acos(dot(vvinfou, vvinfin)/(norm(vvinfin)*norm(vvinfou))); % this should be extremely close to delta_max for TestCase 1, so the error must be vinfAlphaCrank_to_VinfCART function then...


delta = acos(dot(v_inf_minus_enceladus, v_inf_plus_enceladus)/(norm(v_inf_minus_enceladus)*norm(v_inf_plus_enceladus)))


% check the formula that relates the changes in pumb and crack with those
% of deflection angle