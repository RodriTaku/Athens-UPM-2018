function [ Sol, Z, Fval, exitflag, output ] = brachistochrone_optimal_sim( x0, y0, x1, y1, h0, g )
% [ Sol, Z, Fval, exitflag, output ] = brachistochrone_optimal_sim_2( 0, 0, 3, 2, 0.5, 9.8 );

idx = [  1,  5,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 74, 76, 77, 78];
rows = [  2,  3,  4,  5,  6,  9, 10, 11, 12, 13, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27, 30, 31, 32, 33, 34, 37, 38, 39, 40, 41, 44, 45, 46, 47, 48, 51, 52, 53, 54, 55, 58, 59, 60, 61, 62, 65, 66, 67, 68,  1,  2,  3,  4,  5,  3, 10,  4, 11,  1,  5,  6, 10, 11, 12,  1,  6,  9, 10, 11, 12,  1,  2,  7,  1,  2,  8,  1,  2,  9, 10, 17, 11, 18,  1, 12, 13, 17, 18, 19,  1, 13, 16, 17, 18, 19,  1,  6,  7,  9, 14,  1,  6,  8,  9, 15,  1,  6,  9, 16, 17, 24, 18, 25,  1, 19, 20, 24, 25, 26,  1, 20, 23, 24, 25, 26,  1, 13, 14, 16, 21,  1, 13, 15, 16, 22,  1, 13, 16, 23, 24, 31, 25, 32,  1, 26, 27, 31, 32, 33,  1, 27, 30, 31, 32, 33,  1, 20, 21, 23, 28,  1, 20, 22, 23, 29,  1, 20, 23, 30, 31, 38, 32, 39,  1, 33, 34, 38, 39, 40,  1, 34, 37, 38, 39, 40,  1, 27, 28, 30, 35,  1, 27, 29, 30, 36,  1, 27, 30, 37, 38, 45, 39, 46,  1, 40, 41, 45, 46, 47,  1, 41, 44, 45, 46, 47,  1, 34, 35, 37, 42,  1, 34, 36, 37, 43,  1, 34, 37, 44, 45, 52, 46, 53,  1, 47, 48, 52, 53, 54,  1, 48, 51, 52, 53, 54,  1, 41, 42, 44, 49,  1, 41, 43, 44, 50,  1, 41, 44, 51, 52, 59, 53, 60,  1, 54, 55, 59, 60, 61,  1, 55, 58, 59, 60, 61,  1, 48, 49, 51, 56,  1, 48, 50, 51, 57,  1, 48, 51, 58, 59, 66, 60, 67,  1, 61, 62, 66, 67, 68,  1, 62, 65, 66, 67, 68,  1, 55, 56, 58, 63,  1, 55, 57, 58, 64,  1, 55, 58, 65, 68,  1, 62, 63, 65,  1, 62, 64, 65,  1, 62, 65, 69];
cols = [  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  4,  4,  5,  5,  5,  5,  5,  5,  6,  6,  6,  6,  6,  6,  7,  7,  7,  8,  8,  8,  9,  9,  9, 10, 10, 11, 11, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 18, 18, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 25, 25, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 30, 30, 30, 30, 31, 31, 32, 32, 33, 33, 33, 33, 33, 33, 34, 34, 34, 34, 34, 34, 35, 35, 35, 35, 35, 36, 36, 36, 36, 36, 37, 37, 37, 37, 38, 38, 39, 39, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 41, 42, 42, 42, 42, 42, 43, 43, 43, 43, 43, 44, 44, 44, 44, 45, 45, 46, 46, 47, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, 48, 49, 49, 49, 49, 49, 50, 50, 50, 50, 50, 51, 51, 51, 51, 52, 52, 53, 53, 54, 54, 54, 54, 54, 54, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 57, 57, 57, 57, 57, 58, 58, 58, 58, 59, 59, 60, 60, 61, 61, 61, 61, 61, 61, 62, 62, 62, 62, 62, 62, 63, 63, 63, 63, 63, 64, 64, 64, 64, 64, 65, 65, 65, 65, 66, 67, 67, 67, 67, 68, 68, 68, 68, 69, 69, 69, 69];

N = 11;
x = linspace(x0,x1,N).';
y = linspace(y0,y1,N).';
v = linspace(0,sqrt((x1-x0)/(h0*(N-1))^2 + (y1-y0)/(h0*(N-1))^2),N).';
px = ones(N,1);
py = ones(N,1);
pv = ones(N,1);
u = zeros(N,1);
px(1) = nan;
py(1) = nan;
pv(1) = nan;
u(N) = nan;

Sol0 = [x.';y.';v.';u.';px.';py.';pv.'];
Sol0 = [h0;Sol0(:)];
Z0 = Sol0(idx);
options = optimoptions('fsolve','Display','off','Jacobian','on', 'Algorithm', 'levenberg-marquardt','TolX',1e-10,'TolFun',1e-12, 'MaxIter',1000);
[ Z, Fval, exitflag, output ] = brachistochrone_optimal_loop( Z0, x, y, v, g, rows, cols, options );

Sol = Sol0;
Sol(idx) = Z;
Sol = reshape(Sol(2:end),[7,N]);
end