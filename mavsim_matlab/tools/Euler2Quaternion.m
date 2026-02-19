
function quaternion = Euler2Quaternion(phi, theta, psi)
    % converts euler angles to a quaternion
    quaternion = [
                    cos(psi/2)*cos(theta/2)*cos(phi/2) + sin(psi/2)*sin(theta/2)*sin(phi/2);
                    cos(psi/2)*cos(theta/2)*sin(phi/2) - sin(psi/2)*sin(theta/2)*cos(phi/2);
                    cos(psi/2)*sin(theta/2)*cos(phi/2) + sin(psi/2)*cos(theta/2)*sin(phi/2);
                    sin(psi/2)*cos(theta/2)*cos(phi/2) - cos(psi/2)*sin(theta/2)*sin(phi/2)
                ];
end
