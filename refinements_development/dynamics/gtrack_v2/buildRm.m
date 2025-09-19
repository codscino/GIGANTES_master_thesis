function Rm = buildRm(r_enceladus,v_flyby_body)
    % r_encelauds, v_flyby_body row vectors

    b1 = -r_enceladus./norm(r_enceladus);
    b3 = cross(r_enceladus, v_flyby_body)./norm(cross(r_enceladus, v_flyby_body));
    b2 = cross(b3,b1);
    Rm = [ b1' b2' b3' ]';
    
end