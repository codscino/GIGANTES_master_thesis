function Rm = buildRm_robust(r_sc, h_ref)
    % Builds a continuous and robust rotation matrix by using a stable
    % reference vector (h_ref), typically the orbital angular momentum vector
    % of the flyby body, to avoid singularities.

    TOL = 1e-9;

    % b1 is the primary axis, pointing from the spacecraft to Saturn's center.
    % It is always well-defined.
    b1 = -r_sc(:)' / norm(r_sc); % Ensure it's a row vector

    % Define the second basis vector, b2, by taking the cross product of the
    % stable reference vector (h_ref) and our primary axis (b1). This defines
    % a direction in the plane perpendicular to the reference vector.
    b2 = cross(h_ref, b1);

    % Check for the rare edge case where r_sc is aligned with h_ref.
    if norm(b2) < TOL
        % If they are aligned, h_ref can't be used to define the plane.
        % We use a fixed, alternate inertial vector as a fallback.
        ref_alt = [0, 1, 0]; % Inertial Y-axis
        b2 = cross(ref_alt, b1);
    end
    b2 = b2 / norm(b2); % Normalize b2 to be a unit vector

    % The third basis vector, b3, completes the right-handed system.
    % Since b1 and b2 are orthogonal unit vectors, b3 will be one too.
    b3 = cross(b1, b2);

    % Construct the final rotation matrix. The basis vectors are the rows.
    Rm = [b1; b2; b3];
end