%If we have Matlab's symbolic toolbox we can run the following
syms g h

x = sym('x_%d_',[2,1]);
y = sym('y_%d_',[2,1]);
v = sym('v_%d_',[2,1]);
u = sym('u_%d_',[1,1]);

EQS = @(x,y,v,u) ...
[       x(2) - x(1) - h*v(1)*sin(u);
        y(2) - y(1) - h*v(1)*cos(u);
   v(2) - v(1) - h*g*cos(u)];

px = sym('px_%d_',[2,1]);
py = sym('py_%d_',[2,1]);
pv = sym('pv_%d_',[2,1]);
P_vec = [px(2),py(2),pv(2)];

C = @(x,y,v,u) ones(size(u));
J = h*(C(x,y,v,u)) + P_vec*EQS(x,y,v,u);

%The result should be:
%J =  h - pv(2).*(v(1) - v(2) + g*h*cos(u(1))) - py(2).*(y(1) - y(2) + h*v(1).*cos(u(1))) - px(2).*(x(1) - x(2) + h*v(1).*sin(u(1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms g h
N = 11;
x = sym('x_%d_',[N,1]);
y = sym('y_%d_',[N,1]);
v = sym('v_%d_',[N,1]);
px = sym('px_%d_',[N,1]);
py = sym('py_%d_',[N,1]);
pv = sym('pv_%d_',[N,1]);
u = sym('U_%d_%d_',[N,1]);

J = @(x,y,v,u,px,py,pv) h - pv(2).*(v(1) - v(2) + g*h*cos(u(1))) - py(2).*(y(1) - y(2) + h*v(1).*cos(u(1))) - px(2).*(x(1) - x(2) + h*v(1).*sin(u(1)));
S = sym(zeros(N-1,1));
for i = 1:N-1
    S(i) = J(x(i:i+1),y(i:i+1),v(i:i+1),u(i:i+1),px(i:i+1),py(i:i+1),pv(i:i+1));
end
S = sum(S);

var_vec = [px.';py.';pv.';u.';x.';y.';v.'];
var_vec = [h;var_vec(:)];
unused_vec = [px(1);py(1);pv(1);u(N,:).'];
fixed_vec = [x(1);y(1);v(1);x(N);y(N)];
[~, idx] = setdiff(var_vec,union(unused_vec,fixed_vec));
un_vec = var_vec(sort(idx));

var_vec_2 = [x.';y.';v.';u.';px.';py.';pv.'];
var_vec_2 = [h;var_vec_2(:)];
[~, idx_2] = setdiff(var_vec_2,union(unused_vec,fixed_vec));
un_vec_2 = var_vec_2(sort(idx_2));

dS = gradient(S, un_vec);
ddS = jacobian(dS, un_vec_2);

% This all goes inside our brachistochrone simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ Sol, Z, Fval, exitflag, output ] = brachistochrone_optimal_sim( 0, 0, 3, 2, 0.5, 9.8 );