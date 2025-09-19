function fig = plotMoons(IDS, idcentral, holdon)

if nargin == 2
    holdon = 1;
elseif nargin == 3
    if isempty(holdon)
        holdon = 0;
    end
end

% --> gravitational parameter of the central body
[mu_planet, ~] = centralBody_Planet(idcentral);

IDS = unique(IDS);

t = date2mjd2000([ 2023 1 1 0 0 0 ]);

if holdon == 1
    fig = gcf;
    fig.Color = [ 1 1 1 ];
    hold on; grid on; axis equal; 
    xlabel( 'x [km]' ); ylabel( 'y [km]' ); zlabel( 'z [km]' );
else
    fig = figure( 'Color', [1 1 1] );
    hold on; grid on; axis equal; 
    xlabel( 'x [km]' ); ylabel( 'y [km]' ); zlabel( 'z [km]' );
end

[~, namecentral] = planetsMoonsNames(1, idcentral);

plot3( 0, 0, 0, 'o', ...
    'MarkerEdgeColor', 'Black', 'MarkerFaceColor', 'Yellow',...
    'MarkerSize', 10,...
    'DisplayName', namecentral );

for inds = 1:length(IDS)

    [rr, vv, kep] = approxEphem_CC(IDS(inds), t, idcentral);
    Tm            = 2*pi*sqrt( kep(1)^3 / mu_planet );
    [~,yy0]      = propagateKepler_tof(rr, vv, Tm, mu_planet);

    name = planetsMoonsNames(IDS(inds), idcentral);
    plot3( yy0(:,1), yy0(:,2), yy0(:,3), 'LineWidth', 2, 'DisplayName', name );
    % plot3( yy0(:,1), yy0(:,2), yy0(:,3), 'LineWidth', 2, 'HandleVisibility', 'Off' );

end

end