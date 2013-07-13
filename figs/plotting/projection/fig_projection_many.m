function fig_coordinates

A = [0,1,0;0,0,0;0,0,0];
lambda = 5;
B = 1/(1+lambda^2)*(lambda^2 * A - A');
Tjeff = 2*pi*(1+lambda^2)/lambda;
num_times = 200;
ts = linspace(0,Tjeff,num_times);
q05 = [10;0;0];
q01 = [0;1.1;0.4];
q02 = [7;0;1.5];
q03 = [0.5;0;1.5];
q04 = [3;0;1.5];


qs1 = zeros(3,num_times);
ns1 = zeros(3,num_times);

qs2 = zeros(3,num_times);
ns2 = zeros(3,num_times);

qs3 = zeros(3,num_times);
ns3 = zeros(3,num_times);

qs4 = zeros(3,num_times);
ns4 = zeros(3,num_times);

qs5 = zeros(3,num_times);
ns5 = zeros(3,num_times);

for k=1:num_times
   q1 = expm(B*ts(k))*q01;
   n1 = q1/norm(q1);
   qs1(:,k) = q1;
   ns1(:,k) = n1;

   q2 = expm(B*ts(k))*q02;
   n2 = q2/norm(q2);
   qs2(:,k) = q2;
   ns2(:,k) = n2;

   q3 = expm(B*ts(k))*q03;
   n3 = q3/norm(q3);
   qs3(:,k) = q3;
   ns3(:,k) = n3;

   q4 = expm(B*ts(k))*q04;
   n4 = q4/norm(q4);
   qs4(:,k) = q4;
   ns4(:,k) = n4;

   q5 = expm(B*ts(k))*q05;
   n5 = q5/norm(q5);
   qs5(:,k) = q5;
   ns5(:,k) = n5;
end

%plot3(qs(1,:),qs(2,:),qs(3,:)); hold all
%plot3(ns(1,:),ns(2,:),ns(3,:)); 


do_text=0;
fig = figure('position',[100,100,1000,1000]);
ax = subaxis(1,1,1,'margin',0,'padding',0);

colors = [170,39,61; ...
          39,61,170; ...
          ]/255;

hold(ax);

set(ax, 'visible','off');
set(fig,'renderer','opengl');
set(fig,'color',[1 1 1]);

%plane
plane_x = 5*[1, -1, -1, 1];
plane_y = 5*[1, 1, -1, -1];
plane_z = q01(3)*[1, 1, 1,1];
plane_c = [0, 0, 0, 0];
%fill3(plane_x, plane_y, plane_z, plane_c, ...
%    'FaceColor', [.9,.9,.9], ...
%    'FaceAlpha', .6, ...
%    'EdgeColor', [.5,.5,.5], ...
%    'LineWidth', 5 ...
%    );


color = colors(1,:);

%[x,y,z] = tubeplot(qs1,...
%0.01, ... % radius
%24);
%surf(x,y,z,...
%'FaceColor', color,...
%'LineStyle', 'none');    


%sphere
[x,y,z] = sphere(100);
surf(x,y,z,'LineStyle','none', 'facecolor',[.6,.6,.6],'facealpha',.7);

[x,y,z] = tubeplot(ns1,...
0.01, ... % radius
24);
surf(x,y,z,...
'FaceColor', colors(2,:),...
'LineStyle', 'none');    

[x,y,z] = tubeplot(ns2,...
0.01, ... % radius
24);
surf(x,y,z,...
'FaceColor', colors(2,:),...
'LineStyle', 'none');    

[x,y,z] = tubeplot(ns3,...
0.01, ... % radius
24);
surf(x,y,z,...
'FaceColor', colors(2,:),...
'LineStyle', 'none');    


[x,y,z] = tubeplot(ns4,...
0.01, ... % radius
24);
surf(x,y,z,...
'FaceColor', colors(2,:),...
'LineStyle', 'none');    


[x,y,z] = tubeplot(ns5,...
0.01, ... % radius
24);
surf(x,y,z,...
'FaceColor', colors(2,:),...
'LineStyle', 'none');  

[x,y,z] = tubeplot(-ns4,...
0.01, ... % radius
24);
surf(x,y,z,...
'FaceColor', colors(2,:),...
'LineStyle', 'none');  

[x,y,z] = tubeplot(-ns1,...
0.01, ... % radius
24);
surf(x,y,z,...
'FaceColor', colors(2,:),...
'LineStyle', 'none');    

[x,y,z] = tubeplot(-ns2,...
0.01, ... % radius
24);
surf(x,y,z,...
'FaceColor', colors(2,:),...
'LineStyle', 'none');    

[x,y,z] = tubeplot(-ns3,...
0.01, ... % radius
24);
surf(x,y,z,...
'FaceColor', colors(2,:),...
'LineStyle', 'none');    


e_offset = 0;
e_length = 2;
color = [.1,.1,.1];
mArrow3([0,0,0],[0,0,e_length], ...
    'facealpha', 1, ...
    'color', color, ...
    'stemWidth', 0.01,...
    'tipWidth', 0.03); 
mArrow3([0,0,0],[e_length,0,0], ...
    'facealpha', 1, ...
    'color', color, ...
    'stemWidth', 0.01,...
    'tipWidth', 0.03); 

mArrow3([0,0,0],[0,e_length,0], ...
    'facealpha', 1, ...
    'color', color, ...
    'stemWidth', 0.01,...
    'tipWidth', 0.03); 




material dull;

% axis
axis(4*[-1,1,-1,1,-1,1])
axis equal;

view(134,15);
camlight('headlight', 'infinite')
camlight('left', 'infinite')

lighting flat


aafig = myaa()
close(fig);



end

function [x,y,z,c] = get_ellipsoid_meshdata(r, n, lambda)

        particle_size =.5;
        R = vrrotvec2mat(vrrotvec([0,0,1],n));
        rx = particle_size;
        ry = particle_size;
        rz = particle_size;
        if lambda > 1
            rx = rx / lambda;
            ry = ry / lambda;
        else
            rz = rz*lambda;
        end
        [x,y,z] = ellipsoid(0,0,0,rx, ry, rz,100);
        c = double(z > 0) + 1;
        points = [x(:), y(:), z(:)]';
        points = R*points;

        x = r(1) + reshape(points(1,:), size(x));
        y = r(2) + reshape(points(2,:), size(y));
        z = r(3) + reshape(points(3,:), size(z));


        
end