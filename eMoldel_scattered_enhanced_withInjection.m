close all;
clear all;
%set imulation param 
x_size = 10; %unit = m
y_size = 10;
Ts = 10^(-6); %sampling time (unit = s)
plot_num = 20;
iter_num = 500;
grid_numX = 50;
grid_numY = 50; %E-density and temperature maps grid number

%inserted boxes define by interval of x and y axis [x_L x_R y_down y_up]
box_num = 2;
box = [4 6 6 y_size; 4 6 0 4];
%plot the box
figure(2);
for i=1:box_num
    boxL = box(i, 2) - box(i, 1);
    boxH = box(i, 4) - box(i, 3);
    figure(2);
    rectangle('Position', [box(i, 1) box(i, 3) boxL boxH]);
    axis([0 x_size 0 y_size])
    hold on;
end

%e- parameter
m = 9.11 * 10^(-31);
T = 300; %unit = K
kb = 1.38 * 10^(-23); %Boltzmann constant
v_th = sqrt(3*kb*T / m) %thermal velocity

e_num = 1000; %number of electron
tmn = 0.2 * 10^(-12);
p_scatter = 1 - exp(-Ts / tmn);
avg_temp = zeros([1, iter_num]);

%initial position
e_posx = randi([0 min(box(:, 1))], e_num, 1);
e_posy = randi([0 (y_size-1)], e_num, 1); %generate random position of n electrons withing the area
e_posx_old = zeros(e_num, 1);
e_posy_old = zeros(e_num, 1);

%initial speed
std = sqrt(m / (2*pi*kb*T));
rand_v = normrnd(v_th, std, [1 e_num]);
%plot the speed distribution
figure(1);
hist(rand_v);
e_vx = zeros(e_num, 1); 
e_vy = zeros(e_num, 1); 
for i=1:e_num
    low = round((-60 /180) * pi);
    high = round((60 / 180) * pi);
    dir = randi([low high]);
    e_vx(i) = rand_v(i) * cos(dir); %new Vx
    e_vy(i) = rand_v(i) * sin(dir); %new Vy
end

pause(5);

%iteration start
for i=1:iter_num

    if(i >= 5)

        %store the prev position
        e_posx_old = e_posx;
        e_posy_old = e_posy; 
        %position update
        e_posx = e_posx + (Ts*e_vx);
        e_posy = e_posy + (Ts*e_vy);
        
        for j=1:e_num
            vx_old = e_vx(j);
            vy_old = e_vy(j);
            e_vx(j) = vx_old;
            e_vy(j) = vy_old;
    
            %check y-axis
            if (e_posy(j) >= y_size) || (e_posy(j) <= 0)
                e_vy(j) = -vy_old;
            else
                e_vy(j) = vy_old;
            end

            %check x-axis (same as y boundary now
            if (e_posx(j) >= x_size) || (e_posx(j) <= 0)
                e_vx(j) = -vx_old;
            else
                e_vx(j) = vx_old;
            end
           
            %check inserted boxes
            for b=1:box_num
                if (e_posx(j) >= box(b, 1)) && (e_posx(j) <= box(b, 2)) && (e_posy(j) >= box(b, 3)) && (e_posy(j) <= box(b, 4)) %if particl in the dead region
                    if (e_posx_old(j) <= box(b, 1) || (e_posx_old(j) >= box(b, 2))) %if partivle come from the sides of the box
                        e_vx(j) = -vx_old;
                    elseif ((e_posy_old(j) < box(b, 3)) || (e_posy_old(j) > box(b, 4))) %if particle top/bottom of the box
                        e_vy(j) = -vy_old;
                    else
                        e_vx(j) = -vx_old;
                        e_vy(j) = -vy_old;
                    end
                end
            end
    
            %update position
            e_posy(j) = e_posy(j) + (Ts*e_vy(j));
            e_posx(j) = e_posx(j) + (Ts*e_vx(j));

            %{
            %check x boundary
            if (e_posx(j) >= x_size)
                e_posx(j) = 0;
                e_posx_old(j) = 0;
            elseif (e_posx(j) < 0)
                e_posx(j) = x_size;
                e_posx_old(j) = x_size;
            end 
            %}

    
            %speed for next iteration
            if (randi([0 100]) <= p_scatter)
                dir = randi([0 360]);
                e_vx(j) = v_th * cos(dir); %new Vx
                e_vy(j) = v_th * sin(dir); %new Vy
            end
        end
        
        %plotting
        traj_x = [e_posx_old e_posx];
        traj_y = [e_posy_old e_posy];
        for k=1:plot_num
            figure(2);
            plot(traj_x(k, :), traj_y(k, :), 'SeriesIndex', k);
            hold on;
        end
    
        %average temperature calculation
        ev = sqrt((e_vx.^2) + (e_vy.^2));
        ev_avg = sum(ev) / e_num;
        avg_temp(i) = ((ev_avg^2) * m) / (3*kb);
        
        pause(0.0001); %change this number to change the animation speed

    end
end

%plot avg temperature
figure(3);
plot(linspace(0, iter_num, iter_num), avg_temp, ".-")

%plot the E-density and termperature map
grid_x = x_size / grid_numX;
grid_y = y_size / grid_numY;
e_dens = zeros(grid_numY, grid_numX);
t_map = zeros(grid_numY, grid_numX);
for i=0:(grid_numY - 1)
    for j=0:(grid_numX - 1)
        for k=1:e_num
            up_limx = (i*grid_x) + grid_x;
            low_limx = i*grid_x;
            up_limy = (j*grid_y) + grid_y;
            low_limy = j*grid_y;
            if (e_posx(k) < up_limx) && (e_posx(k) >= low_limx) && (e_posy(k) < up_limy) && (e_posy(k) >= low_limy)
                e_dens(i+1, j+1) = e_dens(i+1, j+1) + 1;
                ev = sqrt((e_vx(k)^2) + (e_vy(k)^2));
                t_map(i+1, j+1) = t_map(i+1, j+1) + (((ev^2) * m) / (3*kb));
            end 
        end
    end
end
denX = linspace(0, grid_numX, grid_numX);
denY = linspace(0, grid_numY, grid_numY);
[denx deny] = meshgrid(denX, denY);
%E-density
figure(4);
surf(denx, deny, e_dens);
%T map
figure(5);
surf(denx, deny, t_map);




