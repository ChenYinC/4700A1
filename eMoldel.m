close all;
clear all;
%set area (unit = m)
x_size = 3;
y_size = 3;
plot_num = 20;
iter_num = 200; %iteration time

%e- parameter
m = 9.11 * 10^(-31);
T = 300; %unit = K
kb = 1.38 * 10^(-23); %Boltzmann constant
v_th = sqrt(3*kb*T / m) %thermal velocity
Ts = 10^(-6); %sampling time (unit = s)
t_total = 0.1; %total simulation time (unit = s)
e_num = 5000; %number of electron
avg_temp = zeros([1, iter_num]);

%initial value setup
e_posx = randi([0 (x_size-1)], e_num, 1);
e_posy = randi([0 (y_size-1)], e_num, 1); %generate random position of n electrons withing the area
e_posx_old = zeros(e_num, 1);
e_posy_old = zeros(e_num, 1);

%initial velocity of e
e_vx = zeros(e_num, 1); 
e_vy = zeros(e_num, 1); 
for i=1:e_num
    dirx = randi([0 360]);
    diry = randi([0 360]);
    e_vx(i) = v_th * cos(dirx); 
    e_vy(i) = v_th * sin(diry);
end

%iteration start
for i=1:iter_num
    for j=1:e_num        
        %store the prev position
        e_posx_old(j) = e_posx(j);
        e_posy_old(j) = e_posy(j); 
        %position update
        if (e_posy(j) > y_size) | (e_posy(j) < 0)
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
    ev = v_th;
    ev_avg = v_th;
    avg_temp(i) = ((ev_avg^2) * m) / (3*kb);

    %pause(0.0001); %change this number to change the animation speed
end

%plot the temperature line
figure(2);
plot(linspace(0, iter_num, iter_num), avg_temp, ".-")



