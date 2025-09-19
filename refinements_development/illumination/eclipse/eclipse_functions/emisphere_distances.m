function d_mean = emisphere_distances(d, R)
    % d: Distance between the centers of the spheres
    % R: Radius of sphere A

    % d_mean -> consideirng the emisphere of sphere A facing sphere B,
    % d_mean is the average distance of the misphere point to the center
    % of sphere B

    % the formula has been calculated with a double integral on the
    % emisphere

    
    num = (d^2 + R^2)^(3/2) - (d - R)^3;
    den = 3*d*R;
    
    d_mean =  num/den;
end