function [ Sol, statistics ] = cart_n_pendulum_sim( t0, q0, v0, ffun, ufun, dufun, parameters, h, nsteps, varargin )
%parameters = struct('M', 1/2, 'm', 1/4, 'L', 1, 'g', 3/4);

Jrow = [ 1, 5, 2, 4, 5, 7, 8, 5, 8, 5, 8, 1, 3, 4, 5, 6, 7, 8, 1, 4, 5, 7, 8, 2, 3, 4, 5, 6, 7, 8, 2, 4, 5, 7, 8];
Jcol = [ 1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8];

	nvarargs = length(varargin);
    if nvarargs == 0
        max_iter = 400;
    elseif nvarargs == 1
        max_iter = varargin{1};
    else
        error('Unexpected number of arguments')
    end

%Initialize Newton-Raphson
options = optimoptions('fsolve','Jacobian','on','Display','off','MaxIter',max_iter,'TolX',1e-12,'TolFun',1e-14);
statistics = zeros(3,nsteps+1);

%Initialize variable container
T = t0 + linspace(0,nsteps,nsteps+1)*h;
X = zeros(4,nsteps+1);
F = zeros(2,nsteps+1);
U = zeros(2,nsteps+1);
X(:,1) = [q0;v0];
indsXY = 1:4;
Y = zeros(8,1);
Y(indsXY) = X(:,1);

%Main loop
for j = 2:nsteps+1
    %Newton-Raphson iteration
    tic
    [ Y, G, W, function_value, ~, output ] = cart_n_pendulum_sim_loop( Y, X(:,j-1), T(j-1:j), ffun, ufun, dufun, parameters, Jrow, Jcol, h, options );
    statistics(1,j) = toc;
    statistics(2,j) = norm(function_value, 2);
    statistics(3,j) = output.iterations;
    X(:,j) = Y(indsXY);
	F(:,j) = G;
	U(:,j) = W;
end

Sol = struct('T', T, 'X', X, 'F', F, 'U', U);