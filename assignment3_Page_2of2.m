%Spencer Tigere 101001717

%This section of codes finds the particle tractories for parts 2 c and 3 c
%This code also fid the Density plot for Part 3 a and the Average Current
%vs Bottleneck width in part 2 b). Execute the code twice to observe figures for 3 b
%and 3 c
clc
clear
set(0,'DefaultFigureWindowStyle','docked')

global C
global Ecount
global Vx Vy Vtotal x y
Ecount =1000; 
C.mo = 9.10938215e-31;
C.k = 1.3806504e-23; 
electron_charge = -1.60217662e-19;
Temperature =300;
effective_m = 0.26*C.mo;
Length = 200e-9;
Width = 100e-9; 
Thermal_v = sqrt((2*C.k*Temperature)/effective_m);
time_step = 10e-15; 
frame = 100*time_step; 
x = zeros(Ecount, 2);  
y = zeros(Ecount, 2);  
Temperature = zeros(1,2); 
Time = 0;
VisibleEcount = 50; 
tmn = 0.2e-12;
PScat = 1 - exp(-time_step/tmn);
V_Histogram = zeros(Ecount, 1);
bottleneck_X = [80e-9 80e-9 120e-9 120e-9 80e-9];

bottleneck_Y1 = [100e-9 60e-9 60e-9 100e-9 100e-9];
bottleneck_Y2 = [40e-9 0 0 40e-9 40e-9]; 
Specular = true; 
Inside_Box = true; 
Mapping_S = 10e-9;
Density_Mapping = zeros(Width/Mapping_S, Length/Mapping_S);
Temperature_Mapping = zeros(Width/Mapping_S, Length/Mapping_S);

wid_x = 30;
len_y = 20;
change_x = Length/wid_x;
change_y = Width/len_y;
conduction_outside = 1;
conduction_inside = 01e-2;
conductivity = zeros(wid_x,len_y);
G_matrix = sparse (wid_x*len_y, wid_x*len_y);
V_matrix = zeros(1, wid_x*len_y);
Voltage_x = 0.1;

for i = 1:wid_x
    for j = 1:len_y
        n = j + (i - 1)*len_y;
        nxm = j + ((i-1) - 1)*len_y;
        nxp = j + ((i+1) - 1)*len_y;
        nym = (j-1) + (i - 1)*len_y;
        nyp = (j+1) + (i - 1)*len_y;
        if (i > (0.3*wid_x) || i < (0.6*wid_x)) && (j > (0.6*len_y) || j < (0.3*len_y))
            conductivity(i,j) = conduction_inside;
        else
            conductivity(i,j) = conduction_outside;
        end
        
    end
end
for i = 1:wid_x
    for j = 1:len_y
        n = j + (i - 1)*len_y;
        nxm = j + ((i-1) - 1)*len_y;
        nxp = j + ((i+1) - 1)*len_y;
        nym = (j-1) + (i - 1)*len_y;
        nyp = (j+1) + (i - 1)*len_y;
        if (i == 1)
        n = j + (i - 1)*len_y;
        nxm = j + ((i-1) - 1)*len_y;
        nxp = j + ((i+1) - 1)*len_y;
        nym = (j-1) + (i - 1)*len_y;
        nyp = (j+1) + (i - 1)*len_y;
            V_matrix(n) = Voltage_x;
            G_matrix(n,n) = 1;
        elseif (i == wid_x)
        n = j + (i - 1)*len_y;
        nxm = j + ((i-1) - 1)*len_y;
        nxp = j + ((i+1) - 1)*len_y;
        nym = (j-1) + (i - 1)*len_y;
        nyp = (j+1) + (i - 1)*len_y;
            V_matrix(n) = 0;
            G_matrix(n,n) =1;
        elseif (j == 1)
        n = j + (i - 1)*len_y;
        nxm = j + ((i-1) - 1)*len_y;
        nxp = j + ((i+1) - 1)*len_y;
        nym = (j-1) + (i - 1)*len_y;
        nyp = (j+1) + (i - 1)*len_y;
            G_matrix(n,n) = -((conductivity(i,j) + conductivity(i-1,j))/2) - ((conductivity(i,j) + conductivity(i+1,j))/2) - ((conductivity(i,j) + conductivity(i,j+1))/2);
            G_matrix(n, nxm) = ((conductivity(i,j) + conductivity(i-1,j))/2);
            G_matrix(n,nxp) = ((conductivity(i,j) + conductivity(i+1,j))/2);
            G_matrix(n, nyp) = ((conductivity(i,j) + conductivity(i,j+1))/2);
        elseif (j == len_y)
        n = j + (i - 1)*len_y;
        nxm = j + ((i-1) - 1)*len_y;
        nxp = j + ((i+1) - 1)*len_y;
        nym = (j-1) + (i - 1)*len_y;
        nyp = (j+1) + (i - 1)*len_y;
            G_matrix(n,n) = -((conductivity(i,j) + conductivity(i-1,j))/2) - ((conductivity(i,j) + conductivity(i+1,j))/2) - ((conductivity(i,j) + conductivity(i,j-1))/2);
            G_matrix(n,nxm) = ((conductivity(i,j) + conductivity(i-1,j))/2);
            G_matrix(n,nxp) = ((conductivity(i,j) + conductivity(i+1,j))/2);
            G_matrix(n,nym) = ((conductivity(i,j) + conductivity(i,j-1))/2);
        else
        n = j + (i - 1)*len_y;
        nxm = j + ((i-1) - 1)*len_y;
        nxp = j + ((i+1) - 1)*len_y;
        nym = (j-1) + (i - 1)*len_y;
        nyp = (j+1) + (i - 1)*len_y;
            G_matrix(n,n) = -((conductivity(i,j) + conductivity(i-1,j))/2) - ((conductivity(i,j) + conductivity(i+1,j))/2) - ((conductivity(i,j) + conductivity(i,j-1))/2) - ((conductivity(i,j) + conductivity(i,j+1))/2);
            G_matrix(n,nxm) = ((conductivity(i,j) + conductivity(i-1,j))/2);
            G_matrix(n,nxp) = ((conductivity(i,j) + conductivity(i+1,j))/2);
            G_matrix(n,nym) = ((conductivity(i,j) + conductivity(i,j-1))/2);
            G_matrix(n,nyp) = ((conductivity(i,j) + conductivity(i,j+1))/2);
        end
    end
end

Solution = G_matrix\V_matrix';
surface = zeros(wid_x,len_y);
for i = 1:wid_x
    for j = 1:len_y
        n = j + (i - 1)*len_y;
        nxm = j + ((i-1) - 1) * len_y;
        nxp = j + ((i+1) - 1) * len_y;
        nym = (j-1) + (i - 1) * len_y;
        nyp = (j+1) + (i - 1) * len_y;
        surface(i,j) = Solution(n);
    end
end
[Electricfield_x, Electricfield_y] = gradient(-surface);
Force_x = electron_charge*Electricfield_x;
Force_y = electron_charge*Electricfield_y;
Acceleration_x = Force_x /effective_m;
Acceleration_y = Force_y /effective_m;

for i = 1:Ecount
    x(i,1) = rand()*200e-9;
    y(i,1) = rand()*100e-9;
    Inside_Box = true;
    while Inside_Box == true
        if (x(i) >= 40e-9 && x(i) <= 120e-9) && (y(i) >= 60e-9 ||...
                y(i) <= 40e-9)
            x(i,1) = rand * 200e-9;
            y(i,1) = rand * 100e-9;
        else
            Inside_Box = false;
        end
    end
    
end

for i = 1:Ecount
    
Vx(1:Ecount) = Thermal_v * randn;
Vy(1:Ecount) = Thermal_v * randn;
end

figure(8)
subplot(2,1,1);
plot(bottleneck_X, bottleneck_Y1, bottleneck_X, bottleneck_Y2)
axis([0 Length 0 Width]);
title('Particle Trajecories');
xlabel('x');
ylabel('y');
hold on;

while Time < frame
    subplot(2,1,1)
    for j = 1:Ecount
        leaking = true;
        if PScat> rand
                Vx(j) = Thermal_v * randn;
                Vy(j) = Thermal_v * randn;
        end
        x_index = round((x(j,2)/Length) * 30);
        y_index = round((y(j,2)/Width)*20);
        if x_index < 1
            x_index = 1;
        elseif x_index > 30 
                x_index = 30;
        end
        if y_index < 1
            y_index = 1;
        elseif y_index > 20
            y_index = 20;
        end
        
        Vx(j) =  Vx(j) + Acceleration_x(x_index,y_index)*time_step;
        Vy(j) =  Vy(j) + Acceleration_y(x_index,y_index)*time_step;
        x(j,2) = x(j,1);
        y(j,2) = y(j,1);
        x(j,1) = x(j,1) + (time_step * Vx(j));
        y(j,1) = y(j,1) + (time_step * Vy(j));
        
        if (x(j,1) >= 80e-9 && x(j,1) <= 120e-9) && y(j,1) >= 60e-9
            
                if y(j,2) < 60e-9
                    Vy(j) = -Vy(j);
                    y(j,1) = 60e-9;
                    y(j,2) = 60e-9;
                    
                elseif x(j,2) < 80e-9
                    Vx(j) = -Vx(j);
                    x(j,1) = 80e-9;
                    x(j,2) = 80e-9;
                    
                elseif x(j,2) > 120e-9
                    Vx(j) = -Vx(j);
                    x(j,1) = 120e-9;
                    x(j,2) = 120e-9;
                end
                
            if Specular == true
                
                x(j,1) = x(j,2) + Vx(j)*time_step;
                y(j,1) = y(j,2) + Vy(j)*time_step;
            else
                
             Vx(j) = Thermal_v * randn;
             Vy(j) = Thermal_v * randn;
             
             while leaking == true
                 if(x(j,2) < 80e-9 && Vx(j) >= 0) || ...
                         (x(j,2) > 120e-9 && Vx(j) <= 0) || ...
                         (y(j,2) < 60e-9 && Vy(j) >= 0)      
                     Vx(j) = Thermal_v * randn;
                     Vy(j) = Thermal_v * randn;
                 else
                     leaking = false;
                 end
             end
             x(j,1) = x(j,2) + Vx(j)*time_step;
             y(j,1) = y(j,2) + Vy(j)*time_step;
            end
        end
        if (x(j,1) >= 80e-9 && x(j,1) <= 120e-9) && y(j,1) <= 40e-9
                if y(j,2) > 40e-9
                    Vy(j) = -Vy(j);
                    y(j,1) = 40e-9;
                    y(j,2) = 40e-9; 
                elseif x(j,2) < 80e-9
                    Vx(j) = -Vx(j);
                    x(j,1) = 80e-9;
                    x(j,2) = 80e-9;
                elseif x(j,2) > 120e-9
                    Vx(j) = -Vx(j);
                    x(j,1) = 120e-9;
                    x(j,2) = 120e-9;
                end
            if Specular == true
                x(j,1) = x(j,2) + Vx(j)*time_step;
                y(j,1) = y(j,2) + Vy(j)*time_step;
            else
             Vx(j) = Thermal_v * randn;
             Vy(j) = Thermal_v * randn;
             while leaking == true
                 if(x(j,2) < 80e-9 && Vx(j) >= 0) || ...
                         (x(j,2) > 120e-9 && Vx(j) <= 0) || ...
                         (y(j,2) > 40e-9 && Vy(j) <= 0)
                     Vx(j) = Thermal_v * randn;
                     Vy(j) = Thermal_v * randn;
                 else
                     leaking = false;
                 end
             end
             x(j,1) = x(j,2) + Vx(j)*time_step;
             y(j,1) = y(j,2) + Vy(j)*time_step;
            end
        end
        if x(j,1) > Length
            x(j,2) = 0;
            x(j,1) = time_step * Vx(j);
        end
        if x(j,1) < 0
            x(j,2) = Length;
            x(j,1) = x(j,2) + (time_step * Vx(j));
        end
        if y(j,1) > Width || y(j,1) < 0
            Vy(j) = -Vy(j);
        end
        XPlot = [x(j,2) x(j,1)];
        YPlot = [y(j,2) y(j,1)];
        if j < VisibleEcount
        plot(XPlot,YPlot);
        end
        
       VTotal = sqrt(Vx(j)^2 + Vy(j)^2);    
    end

    AvgTemperature = Temperature(1,2)/Ecount;
    TemperaturePlot = [Temperature(1,1) AvgTemperature];
    TimePlot = [(Time - time_step) Time];
    subplot(2,1,2);
    plot(TimePlot, TemperaturePlot);
    Temperature(1,1) = AvgTemperature;
    AvgTemperature = 0;
    Temperature(1,2) = 0;
    pause(1e-19)
    Time = Time + time_step;
end 
for i = 1:(Length/Mapping_S)
    for j = 1:(Width/Mapping_S)
        for m = 1:Ecount
            if(x(m,1) > Mapping_S*(i -1)) && ...
                    (x(m,1) < Mapping_S*(i)) && ...
                    (y(m,1) > Mapping_S*(j - 1)) && ...
                    (y(m,1) < Mapping_S*(j))
             
                Vtotal(m) = sqrt(Vx(m)^2 + Vy(m)^2);
                
                Density_Mapping(j, i) = Density_Mapping(j, i) + 1;
                Temperature_Mapping(j, i) = Temperature_Mapping(j,i) + ...
                    (effective_m*Vtotal(m)^2)/(2*C.k);
            end
            Temperature_Mapping(j,i) = Temperature_Mapping(j,i)/Density_Mapping(j,i);
        end
    end
end
figure(9)
imagesc(Density_Mapping)
xlabel('x');
ylabel('y');
set(gca, 'Ydir', 'Normal')
title('Denisty Plot')