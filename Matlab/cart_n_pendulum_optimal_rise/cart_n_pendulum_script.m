N = 40;
X0 = 0; VX0 = 0;
Theta0 = 0; VTheta0 = 0;
X1 = 0; VX1 = 0;
Theta1 = pi; VTheta1 = 0;
g = -0.75; h = 0.3;
M = 0.5; m = 0.25;
l = 1;
sigma = [1,5,10,50,100,100,500,1000];
ulim = 1;

x = linspace(X0,X1,N).';
temp = tanh(linspace(-1,1,N))-tanh(-1);
theta = (temp/temp(end)*(Theta1-Theta0)+Theta0).';
vx = linspace(VX0,VX1,N).';
vtheta = linspace(VTheta0,VTheta1,N).';
px = ones(N,1);
ptheta = ones(N,1);
pvx = zeros(N,1);
pvtheta = zeros(N,1);
X = zeros(N,2);
Theta = zeros(N,2);
VX = zeros(N,2);
VTheta = zeros(N,2);
PX = zeros(N,2);
PTheta = zeros(N,2);
PVX = zeros(N,2);
PVTheta = zeros(N,2);
U = zeros(N,2);
for i = 1:N-1
    X(i,:) = linspace(x(i),x(i+1),2);
    Theta(i,:) = linspace(theta(i),theta(i+1),2);
    VX(i,:) = linspace(vx(i),vx(i+1),2);
    VTheta(i,:) = linspace(vtheta(i),vtheta(i+1),2);
    PVX(i,:) = linspace(pvx(i),pvx(i+1),2);
    PVTheta(i,:) = linspace(pvtheta(i),pvtheta(i+1),2);
end

idx = [   5,   6,   7,   8,   9,  10,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 721, 722, 723, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 743, 744, 745, 746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 765, 766, 767, 768, 769, 770, 771, 772, 773, 774, 775, 776, 777, 778, 779, 780, 781, 782, 787, 788, 789, 790, 791, 792, 793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 821, 822, 823, 824, 825, 826, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 853, 854, 855, 856, 857, 858, 869, 870, 875, 876];
Sol0 = [x.';theta.';vx.';vtheta.';VX.';VTheta.';U.';px.';ptheta.';PX.';PTheta.';pvx.';pvtheta.';PVX.';PVTheta.'];
Z0 = Sol0(:);
Y = Z0(idx);

options = optimoptions('fsolve','Display','off','Jacobian','on','TolX',1e-10,'TolFun',1e-12, 'MaxIter',1000);
for i = 1:8
    [ Y, Fval, exitflag, output ] = cart_n_pendulum( Y, x, theta, vx, vtheta, M, m, l, g, -sigma(i), ulim, h, options );
end

Sol = Sol0;
Sol(idx) = Y;
x = Sol(1,:).';
theta = Sol(2,:).';
vx = Sol(3,:).';
vtheta = Sol(4,:).';
VX = Sol(5:6,:).';
VTheta = Sol(7:8,:).';
U = Sol(9:10,:).';
px = Sol(11,:).';
ptheta = Sol(12,:).';
PX = Sol(13:14,:).';
PTheta = Sol(15:16,:).';
pvx = Sol(17,:).';
pvtheta = Sol(18,:).';
PVX = Sol(19:20,:).';
PVTheta = Sol(21:22,:).';
t = h*linspace(0,N-1,N).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pp_x = spline(t,[vx(1),x.',vx(end)]);
pp_theta = spline(t,[vtheta(1),theta.',vtheta(end)]);
pp_u = spline(t,[0,U(:,1).',0]);
N_new = 300;
t_new = linspace(t(1),t(end),N_new);
x_new = ppval(pp_x, t_new);
theta_new = ppval(pp_theta, t_new);
u_new = ppval(pp_u, t_new);

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
        p_rod = plot([x_new(i), x_new(i)+l*sin(theta_new(i))],[0, -l*cos(theta_new(i))],'g', 'LineWidth',1);
        hold on
        p_control = quiver(x_new(i),0,u_new(i),0,'r');
        hold on
        p_bob = scatter(x_new(i)+l*sin(theta_new(i)),-l*cos(theta_new(i)),25,'MarkerEdgeColor',[0 .5 .5],...
                'MarkerFaceColor',[0 .7 .7],...
                'LineWidth',1.5);
        hold off
        axis([min(x)-2*max(l,box_length/2),max(x)+2*max(l,box_length/2),min(x)-2*max(l,box_length/2),max(x)+2*max(l,box_length/2)])
        daspect([1,1,1])
%         frame = getframe(F); 
%         im = frame2im(frame); 
%         [imind,cm] = rgb2ind(im,256); 
%         % Write to the GIF File 
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 
    else
        set(p_box,'XData',x_box,'YData',y_box)
        set(p_rod,'XData',[x_new(i), x_new(i)+l*sin(theta_new(i))],'YData',[0, -l*cos(theta_new(i))])
        set(p_control,'XData',x_new(i),'UData',u_new(i))
        set(p_bob,'XData',x_new(i)+l*sin(theta_new(i)),'YData',-l*cos(theta_new(i)))
        drawnow
        frame = getframe(F); 
%         im = frame2im(frame); 
%         [imind,cm] = rgb2ind(im,256); 
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
    end
    %pause(0.01)
end
