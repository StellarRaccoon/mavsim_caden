
function [phi, theta, psi] = Quaternion2Euler(quaternion)
    % converts a quaternion attitude to an euler angle attitude
    e0 = quaternion(1);
    e1 = quaternion(2);
    e2 = quaternion(3);
    e3 = quaternion(4);
    
    phi = atan2(2*(e0*ex + ey*ez),(e0^2 + ez^2 - ex^2 - ey^2));
    theta = asin(2*(e0*ey - ex*ez));
    psi = atan2(2*(e0*ez + ex*ey),(e0^2 + ex^2 - ey^2 - ez^2));
end
