function [Fig] = Plot_PetalRotation_Trajectory(PetalRotation_Sequence, PetalRotation_FlybyData, pars)


colors = cool(size(PetalRotation_Sequence,1));

% Plot the figure
Fig = figure('Color', [1 1 1]);clf;set(Fig,'defaulttextinterpreter','latex') ;
hold on; grid on; axis equal; 
plotMoons( pars.INPUTS.idMoon, pars.INPUTS.idCentral);
axis equal; hold on;
for i = 1:size(PetalRotation_Sequence, 1)

    PseudoRes = PetalRotation_Sequence(i,:);
    TOF(i)    = PseudoRes(6);
    if i == 1
        [~, rr1, vv1, ~] = vinfAlphaCrank_to_VinfCART(PseudoRes(3), PseudoRes(4), PseudoRes(5) , 0, pars.INPUTS.idMoon, pars.INPUTS.idCentral);
        [~, yy]                  = propagateKepler_tof(rr1, vv1, TOF(i), pars.Planet.mu);
    else
        [~, rr1, vv1, ~] = vinfAlphaCrank_to_VinfCART(PseudoRes(3), PseudoRes(4), PseudoRes(5), sum(TOF(1:i-1))/86400, pars.INPUTS.idMoon, pars.INPUTS.idCentral);
        [~, yy]          = propagateKepler_tof(rr1, vv1, TOF(i), pars.Planet.mu);
    end
    plot3( yy(:,1), yy(:,2), yy(:,3), 'LineWidth', 2,  'Color', colors(i,:));
    plot3(yy(1,1), yy(1,2), yy(1,3), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'MarkerSize', 5 );
    plot3(yy(end,1), yy(end,2), yy(end,3), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green', 'MarkerSize', 5 );

    % Convert last state vector to OE & Store
    OE_SC_end(i,:)  = car2kep(yy(end,:),pars.Planet.mu);
    clear yy;
end

% Retrieve SC True Longitude at the initial flyby
TrueLong_0 = mod(sum(OE_SC_end(1, 4:6)), 2*pi);

% Compute change in True Longitude achieved by each transfer
for i = 2:size(OE_SC_end,1)
        TrueLong(i)         = mod(sum(OE_SC_end(i,4:6)), 2*pi); % [rad]
        TrueLong_Diff(i, :) = TrueLong(i) - TrueLong_0;         % [rad]
end

end
