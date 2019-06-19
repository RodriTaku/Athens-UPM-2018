function [ A, B ] = cart_n_pendulum_linearization( X, param )
	x = X(1);
	theta = X(2);
	vx = X(3);
	vtheta = X(4);

	M = param.M;
	m = param.m;
	g = param.g;
	L = param.L;

	A = [[ 0,                                                                                                                                                                                                                                                      0, 1,                                                     0];
		 [ 0,                                                                                                                                                                                                                                                      0, 0,                                                     1];
		 [ 0,                                                           (m*(2*g*cos(theta)^2 - g + L*vtheta^2*cos(theta)))/(M + m - m*cos(theta)^2) - (m^2*cos(theta)*sin(theta)*(2*L*sin(theta)*vtheta^2 + 2*g*cos(theta)*sin(theta)))/(- m*cos(theta)^2 + M + m)^2, 0,        (2*L*m*vtheta*sin(theta))/(m*sin(theta)^2 + M)];
		 [ 0, (2*m*cos(theta)*sin(theta)*(L*m*cos(theta)*sin(theta)*vtheta^2 + g*m*sin(theta) + M*g*sin(theta)))/(L*(- m*cos(theta)^2 + M + m)^2) - (2*(L*m*(2*cos(theta)^2 - 1)*vtheta^2 + g*m*cos(theta) + M*g*cos(theta)))/(L*(2*M + m - m*(2*cos(theta)^2 - 1))), 0, -(2*m*vtheta*sin(2*theta))/(2*M + m - m*cos(2*theta))]];
	 
	B = [                                        0;
												 0;
						1/(M + m - m*cos(theta)^2);
		  -cos(theta)/(L*(M + m - m*cos(theta)^2))];
end

