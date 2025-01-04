function [Irot,Qrot,theta0] = rotate_IQ(I,Q,theta_step_size)
% theta_step_size = 5, 10, 15, 20 
alpha = atan2d(Q',I');
% look where the density of points is largest
Niqangle_pos = nan(180/theta_step_size,1);
for n=1:180/theta_step_size
  Niqangle_pos(n) = size(alpha(alpha>=(n-1)*theta_step_size & alpha<=n*theta_step_size),1);
end

Niqangle_neg = nan(180/theta_step_size,1);
for n=1:180/theta_step_size
  Niqangle_neg(n) = size(alpha(alpha<=-(n-1)*theta_step_size & alpha>=-n*theta_step_size),1);
end

if max(Niqangle_pos)>=max(Niqangle_neg)
    theta_array = (find(Niqangle_pos==max(Niqangle_pos))-1) * theta_step_size;
    else
    theta_array = -(find(Niqangle_neg==max(Niqangle_neg))-1) * theta_step_size;
end

theta0 = theta_array(1);

theta = pi/4 - theta0 * pi/180;

Irot = I*cos(theta) - Q*sin(theta);
Qrot = I*sin(theta) + Q*cos(theta);

end