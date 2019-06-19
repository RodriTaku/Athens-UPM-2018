function [ Y, F, U, Fval, exitflag, output ] = cart_n_pendulum_sim_loop(Y, X, T, Ffun, Ufun, dUfun, param, Jrow, Jcol, h, options)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    [ Y, Fval, exitflag, output ] = fsolve(@my_fun, Y, options);
	F = [Ffun(T(1)), Ffun(T(2))];
	U = [Ufun(T(1), X(1), X(2), Y(5), Y(7)), Ufun(T(2), Y(1), Y(2), Y(6), Y(8))];

    function [ Fun, dFun ] = my_fun( Y )
	
	F = [Ffun(T(1)), Ffun(T(2))];
	U = [Ufun(T(1), X(1), X(2), Y(5), Y(7)), Ufun(T(2), Y(1), Y(2), Y(6), Y(8))];
	dU = [dUfun(T(1), X(1), X(2), Y(5), Y(7)), dUfun(T(2), Y(1), Y(2), Y(6), Y(8))];
    
	Fun = ...
	[                                                                                                                                                                                                                                                Y(1) - X(1) - (Y(5)*h)/2 - (Y(6)*h)/2;
																																																													 Y(2) - X(2) - (Y(7)*h)/2 - (Y(8)*h)/2;
																																			   Y(5)*(param.M + param.m) - X(3)*(param.M + param.m) - (h*(F(1) + U(1)))/2 - X(4)*param.L*param.m*cos(X(2)) + Y(7)*param.L*param.m*cos(X(2));
																																			   Y(6)*(param.M + param.m) - X(3)*(param.M + param.m) - (h*(F(1) + U(1)))/2 - X(4)*param.L*param.m*cos(X(2)) + Y(8)*param.L*param.m*cos(Y(2));
																											 2*Y(3)*(param.M/2 + param.m/2) - 2*X(3)*(param.M/2 + param.m/2) - (h*(F(1) + U(1)))/2 - (h*(F(2) + U(2)))/2 - X(4)*param.L*param.m*cos(X(2)) + Y(4)*param.L*param.m*cos(Y(2));
																			  (h*param.L*(param.g*param.m*sin(X(2)) - cos(X(2))*F(1) + Y(5)*Y(7)*param.m*sin(X(2))))/2 - X(4)*param.L^2*param.m + Y(7)*param.L^2*param.m - X(3)*param.L*param.m*cos(X(2)) + Y(5)*param.L*param.m*cos(X(2));
																			  (h*param.L*(param.g*param.m*sin(X(2)) - cos(X(2))*F(1) + Y(5)*Y(7)*param.m*sin(X(2))))/2 - X(4)*param.L^2*param.m + Y(8)*param.L^2*param.m - X(3)*param.L*param.m*cos(X(2)) + Y(6)*param.L*param.m*cos(Y(2));
	 (param.L*(2*Y(4)*param.L*param.m - h*cos(Y(2))*F(2) - 2*X(4)*param.L*param.m - h*cos(X(2))*F(1) - 2*X(3)*param.m*cos(X(2)) + 2*Y(3)*param.m*cos(Y(2)) + param.g*h*param.m*sin(X(2)) + param.g*h*param.m*sin(Y(2)) + Y(5)*Y(7)*h*param.m*sin(X(2)) + Y(6)*Y(8)*h*param.m*sin(Y(2))))/2];

    dFun = sparse(Jrow, Jcol, [                                                                                                                    1;
																																	  -(h*dU(2,2))/2;
																																				   1;
																													 -Y(8)*param.L*param.m*sin(Y(2));
																									- (h*dU(3,2))/2 - Y(4)*param.L*param.m*sin(Y(2));
																													 -Y(6)*param.L*param.m*sin(Y(2));
							 (param.L*(h*F(2)*sin(Y(2)) - 2*Y(3)*param.m*sin(Y(2)) + param.g*h*param.m*cos(Y(2)) + Y(6)*Y(8)*h*param.m*cos(Y(2))))/2;
																																   param.M + param.m;
																														   param.L*param.m*cos(Y(2));
																														   param.L*param.m*cos(Y(2));
																																   param.L^2*param.m;
																																				-h/2;
																												   param.M + param.m - (h*dU(4,1))/2;
																																	  -(h*dU(4,1))/2;
																																	  -(h*dU(4,1))/2;
																					param.L*param.m*cos(X(2)) + (Y(7)*h*param.L*param.m*sin(X(2)))/2;
																												(Y(7)*h*param.L*param.m*sin(X(2)))/2;
																												(Y(7)*h*param.L*param.m*sin(X(2)))/2;
																																				-h/2;
																																   param.M + param.m;
																																	  -(h*dU(4,2))/2;
																														   param.L*param.m*cos(Y(2));
																												(Y(8)*h*param.L*param.m*sin(Y(2)))/2;
																																				-h/2;
																										   param.L*param.m*cos(X(2)) - (h*dU(5,1))/2;
																																	  -(h*dU(5,1))/2;
																																	  -(h*dU(5,1))/2;
																							param.m*param.L^2 + (Y(5)*h*param.m*sin(X(2))*param.L)/2;
																												(Y(5)*h*param.L*param.m*sin(X(2)))/2;
																												(Y(5)*h*param.L*param.m*sin(X(2)))/2;
																																				-h/2;
																														   param.L*param.m*cos(Y(2));
																																	  -(h*dU(5,2))/2;
																																   param.L^2*param.m;
																												(Y(6)*h*param.L*param.m*sin(Y(2)))/2]);
    end

end