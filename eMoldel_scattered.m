close all;
clear all;
%set simulation param 
x_size = 10; %unit = m
y_size = 10;
Ts = 10^(-6); %sampling time (unit = s)
plot_num = 8; %number of partivle that plot 
total_num = 5000; %total particle number
iter_num = 600;
total_time = iter_num * Ts; %total simulation time 
grid_numX = 50;
grid_numY = 50; %E-density and temperature maps grid number

%e- parameter
m = 9.11 * 10^(-31);
T = 300; %unit = K
kb = 1.38 * 10^(-23); %Boltzmann constant
v_th = sqrt(3*kb*T / m) %thermal velocity

e_num = total_num; %number of electron
tmn = 0.2 * 10^(-12);
p_scatter = 1 - exp(-Ts / tmn);
avg_temp = zeros([1, iter_num]);
scatter_num = 0; %total number of scatter

%initial position
e_posx = randi([0 x_size], e_num, 1);
e_posy = randi([0 y_size], e_num, 1); %generate random position of n electrons withing the area
e_posx_old = zeros(e_num, 1);
e_posy_old = zeros(e_num, 1);

%initial speed
std = sqrt(m / (2*pi*kb*T));
rand_v = normrnd(v_th, std, [1 e_num]);
figure(1);
hist(rand_v);
e_vx = zeros(e_num, 1); 
e_vy = zeros(e_num, 1); 
for i=1:e_num
    dir = randi([0 360]);
    e_vx(i) = rand_v(i) * cos(dir); %new Vx
    e_vy(i) = rand_v(i) * sin(dir); %new Vy
end

%iteration start
for i=1:iter_num
    %store the prev position
    e_posx_old = e_posx;
    e_posy_old = e_posy; 
    %position update
    e_posx = e_posx + (Ts*e_vx);
    e_posy = e_posy + (Ts*e_vy);
    
    for j=1:e_num
        %check y-axis
        if (e_posy(j) >= y_size) | (e_posy(j) <= 0)
            e_vy(j) = -e_vy(j);
            e_posy(j) = e_posy(j) + (Ts*e_vy(j));
            e_posx(j) = e_posx(j) + (Ts*e_vx(j));
        else
            e_posx(j) = e_posx(j) + (Ts*e_vx(j));
            e_posy(j) = e_posy(j) + (Ts*e_vy(j));
        end

        %check x-axis
        if (e_posx(j) >= x_size)
            e_posx(j) = 0;
            e_posx_old(j) = 0;
        elseif (e_posx(j) < 0)
            e_posx(j) = x_size;
            e_posx_old(j) = x_size;
        end 

        %speed for next iteration
        if (randi([0 100]) <= p_scatter)
            scatter_num = scatter_num + 1;
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

%plot the temperature line
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

%MFP calculation
avgV = sum(sqrt((e_vx.^2) + (e_vy.^2))) / e_num; %average velocity at the end
colliFreq = scatter_num / total_time; %collision frequency
scatter_num
tmn = 1 / colliFreq
mfp = avgV / colliFreq

