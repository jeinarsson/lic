function fig_coordinate

    num_frames = 10;
    steps_per_frame = 7;
    trail = steps_per_frame * 30;
    start_frame = trail + 340;
    
    name = 'triax';
    %name = 'symm';
    
    trajectory = load(sprintf('anton/%s.mat',name),'n_x','n_y','n_z');
    
    
    fig = figure('position',[100,100,100,100]);

    
    ax = subaxis(1,1,1,'Holdaxis',1,'Margin',.02);
    %Other items

    set(ax, 'visible','off');
    set(fig,'renderer','opengl');
    set(fig,'color',[1 1 1]);    
    
    lambda = 7;


    r = [0,0,0];
    n =[1,1,0];
    [x,y,z,c] = get_ellipsoid_meshdata(r,n,lambda);
    
    
    
    particle_mesh = surf(x,y,z,c,'LineStyle','none');
    


    hold on;

    trail_idx = (start_frame-trail):start_frame;
    trail_ns = [trajectory.n_x(trail_idx), trajectory.n_y(trail_idx), trajectory.n_z(trail_idx)]';
    
    [trail_mesh_x ,trail_mesh_y ,trail_mesh_z ]= tubeplot(trail_ns,...
    0.01, ... % radius
    24);

    size(trail_mesh_x)
    trail_length = size(trail_mesh_x, 2)
    cs = linspace(0,1,trail_length);
    trail_alpha = repmat(linspace(0,1,trail_length), 25,1);
    size(trail_alpha)
    trail_mesh = surf(trail_mesh_x,trail_mesh_y,trail_mesh_z,'FaceColor',[0.3,0.3,0.3],'AlphaData',trail_alpha,'FaceAlpha','interp','LineStyle', 'none');     
    



%    mArrow3([0,0,0],[0,0,1], ...
%    'facealpha', 1, ...
%    'color', [0,0,0], ...
%    'stemWidth', 0.01,...
%    'tipWidth', 0.03); 
    
    material dull;

    % axis
    axis(1*[-1,1,-1,1,-1,1])
    axis off
    view(36,6);
    camlight('headlight', 'infinite')
    camlight('left', 'infinite')
    camproj('perspective')
    %camproj('orthographic')
    
    lighting flat
    
    %aafig = myaa()
    %close(fig);
    
    
    k = start_frame;
    for frame_idx = 1:num_frames
        
        k = k + steps_per_frame;

        n = [trajectory.n_x(k); trajectory.n_y(k); trajectory.n_z(k)];
        [x,y,z] = get_ellipsoid_meshdata(r,n,lambda);
        set(particle_mesh, 'xdata',x,'ydata',y,'zdata',z);

        trail_idx = (k-trail):k;
        trail_ns = [trajectory.n_x(trail_idx), trajectory.n_y(trail_idx), trajectory.n_z(trail_idx)]';
        [x ,y ,z ]= tubeplot(trail_ns,...
        0.01, ... % radius
        24);        
        set(trail_mesh, 'xdata',x,'ydata',y,'zdata',z);
        
        aafig = myaa();
        
        movefile('myaa_temp_screendump.png', sprintf('output/%s/%d.png',name,frame_idx-1));
        
        close(aafig);
         %pause(0.05)
       
    end
    

    
   close(fig);



end

function [x,y,z,c] = get_ellipsoid_meshdata(r, n, lambda)

        particle_size = 1;
        
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
        
        colors= [
        170,39,61;
        229,142,26;
        ]/255;
        
        [x,y,z] = ellipsoid(0,0,0,rx, ry, rz,100);
        c = double(x > 0);
        cs = zeros([size(z) 3]);
        for k = 1:3
            cs(:,:,k) = colors(1,k)*c + (1-c)*colors(2,k);
        end
        c = cs;
        points = [x(:), y(:), z(:)]';
        points = R*points;

        x = r(1) + reshape(points(1,:), size(x));
        y = r(2) + reshape(points(2,:), size(y));
        z = r(3) + reshape(points(3,:), size(z));


        
end