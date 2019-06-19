t = Sol.T;
x = Sol.X(1,:);
theta = Sol.X(2,:);
vx = Sol.X(3,:);
vtheta = Sol.X(4,:);
f = [Sol.F(1,2),(Sol.F(1,3:end)+Sol.F(2,2:end-1))/2,Sol.F(2,end)];
u = [Sol.U(1,2),(Sol.U(1,3:end)+Sol.U(2,2:end-1))/2,Sol.U(2,end)];

pp_x = spline(t,[vx(1),x,vx(end)]);
pp_theta = spline(t,[vtheta(1),theta,vtheta(end)]);
pp_f = spline(t,f);
pp_u = spline(t,u);

N_new = 500;
t_new = linspace(t(1),t(end),N_new);
x_new = ppval(pp_x, t_new);
theta_new = ppval(pp_theta, t_new);
f_new = ppval(pp_f, t_new);
u_new = ppval(pp_u, t_new);

f_scaling = 2;
u_scaling = 2;
box_length = 0.3;
box_height = 0.2;
F = figure();
%filename = 'cart_n_pendulum_optimal.gif';
for i = 1:N_new
    x_box = [x_new(i)+box_length/2, x_new(i)+box_length/2, x_new(i)-box_length/2, x_new(i)-box_length/2, x_new(i)+box_length/2];
    y_box = [-box_height/2, box_height/2, box_height/2, -box_height/2, -box_height/2];
    if i == 1
        p_box = plot(x_box,y_box,'b', 'LineWidth',1);
        hold on
        p_rod = plot([x_new(i), x_new(i)+param.L*sin(theta_new(i))],[0, -param.L*cos(theta_new(i))],'g', 'LineWidth',1);
        hold on
        p_control = quiver(x_new(i),0,u_new(i)*u_scaling,0,'r');
        hold on
        p_force = quiver(x_new(i)+param.L*sin(theta_new(i)),-param.L*cos(theta_new(i)),f_new(i)*f_scaling,0,'m');
        hold on
        p_bob = scatter(x_new(i)+param.L*sin(theta_new(i)),-param.L*cos(theta_new(i)),25,'MarkerEdgeColor',[0 .5 .5],...
        'MarkerFaceColor',[0 .7 .7],...
        'LineWidth',1.5);
        hold off
        axis([min(x)-2*max(param.L,box_length/2),max(x)+2*max(param.L,box_length/2),min(x)-2*max(param.L,box_length/2),max(x)+2*max(param.L,box_length/2)])
        daspect([1,1,1])
        %frame = getframe(F);
        %im = frame2im(frame);
        %[imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        %imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0);
    else
        set(p_box,'XData',x_box,'YData',y_box)
        set(p_rod,'XData',[x_new(i), x_new(i)+param.L*sin(theta_new(i))],'YData',[0, -param.L*cos(theta_new(i))])
        set(p_control,'XData',x_new(i),'UData',u_new(i)*u_scaling)
        set(p_force,'XData',x_new(i)+param.L*sin(theta_new(i)),'YData',-param.L*cos(theta_new(i)),'UData',f_new(i)*f_scaling)
        set(p_bob,'XData',x_new(i)+param.L*sin(theta_new(i)),'YData',-param.L*cos(theta_new(i)))
        drawnow
        %frame = getframe(F);
        %im = frame2im(frame);
        %[imind,cm] = rgb2ind(im,256);
        %imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);
    end
    %pause(0.01)
end