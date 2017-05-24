function binarydistillation(~)
clc
%This function uses several different parts and ultimately outputs a 
%McCabe Thiele diagram for a BINARY mixture and a plot showing the effect
%of the reflux ratio on the number of stages needed. Though a lot of edits
%have been made, the framework was built off of code by ashishkshwh
%(GitHub username).

%--------------------------------------------------------------------------

%User inputted information

z_F = .4; %z_F is the feed mole fraction
x_D = .7; %x_D is the distillate mole fraction
x_B = .01; %x_B is the bottoms mole fraction
eta = .7; %eta is the stage efficiency
R_by_Rm = 2; %R_by_Rm is the value of R/Rm
l = 0.3; %Liquid fraction in feed

%Data from Perry's Handbook
x = [.05;.1;.2;.3;.4;.5;.6;.7;.8];
y = [.102;.186;.322;.428;.513;.586;.656;.725;.8];

%--------------------------------------------------------------------------
%Calculates the minimum reflux and feed line required for
%distillation of a given feed using McCabe Thiele Method

Y = (y.*(1-x))./((1-y).*x);

%To plot equilbrium curve

x_A = 0:0.05:1;
y_A = yeq(x_A,x,Y);

v = 1 - l;
m_q = -l/v;
c_q = z_F*(1 - m_q);
x_0 = 0.25; %Initial guess to solve q-line and eq curve
q_fun = @(x_var)(m_q*x_var + c_q - yeq(x_var,x,Y));
options = optimoptions('fsolve', 'Display', 'off');
x_ep = fsolve(q_fun, x_0, options); %End point of the q-line
y_ep = m_q*x_ep + c_q; %End point of the q-line

%Minimum reflux from the pinch point
Rm = (x_D - y_ep)/(y_ep - x_ep);
fprintf('The minimum reflux required is %.2f. \n', Rm);

%--------------------------------------------------------------------------

%Calculates the minimum number of stages required for
%distillation of a given feed using McCabe Thiele Method

stages = 0;
x_e = x_D;
stage_plot = [];
while (x_e > x_B)
    stage_plot = [stage_plot; [x_e, x_e]];
    y_temp = stage_plot(end, 2);
    x_temp = xeq(y_temp,y,Y);
    stages = stages + 1;
    stage_plot = [stage_plot; [x_temp, y_temp]];
    x_e = x_temp;
end

    stages = stages - (x_e - x_B)/(x_e - stage_plot(end-1, 1));
    fprintf('Minimum number of stages is %.2f. \n', stages);
%--------------------------------------------------------------------------

%Calculates the effect of reflux ratio on the number of stages

Ref = [1.1:0.15:2];
Ref(:);
s = size(Ref);
for i = 1:s(2)
    R_by_Rm = Ref(i);    
    [stages] = stg(R_by_Rm,Rm,x_D,x_B,c_q,m_q,y,Y,eta);
    stage_arr(i) = stages;
end

%--------------------------------------------------------------------------

%Calculates the number of stages required with a given R/Rm ratio for
%distillation of a given feed using McCabe Thiele Method

[stages,x_rec,y_rec,x_str,y_str,stage_plot] = stg(R_by_Rm,Rm,x_D,x_B,c_q,m_q,y,Y,eta);
fprintf('Actual number of stages is %.0f with reflux and efficiencies taken into account.\n', stages);

%--------------------------------------------------------------------------

%Plotting

figure('Name','Distillation using McCabe Thiele Method')
subplot(1,2,1)
hold on
grid on
axis equal
axis([0 1 0 1])
plot(x_A, y_A); %Plots the equilibrium line
plot([0 1],[0 1], 'k-');%x=y line
xlabel('x_A');
ylabel('y_A');
plot([x_D, x_D], [0, x_D], '--'); %Shows where distillate composition is
plot([x_B, x_B], [0, x_B], '--'); %Shows where bottoms composition is
plot([z_F, z_F, x_ep], [0, z_F, y_ep], '--'); %Shows where feed is
plot([x_rec x_D], [y_rec x_D], 'r'); %Plots the rectifying operating line
plot([x_B x_str], [x_B y_str], 'r') %Plots the stripping operating line
plot(stage_plot(:, 1), stage_plot(:, 2)); %Plots the stages
title('McCabe Thiele Diagram');
subplot(1,2,2)
plot(Ref, stage_arr);
xlabel('R/R_m');
ylabel('Number of stages');
grid on
title('Reflux Effect on Stages');

end

%Function to calculate y_eq for a given x
function y_A = yeq(x_A,x,Y)
    p = polyfit(x,Y,5);
    alpha = polyval(p,x_A);
    y_A = (alpha.*x_A)./(1 - x_A + alpha.*x_A);
end

%Function to calculate x_eq for a given y
function x_A = xeq(y_A,y,Y)
    p = polyfit(y,Y,5);
    alpha = polyval(p,y_A);
    x_A = y_A./(alpha - alpha.*y_A + y_A);
end

%Function to calculate the number of stages for a given R/Rm
function [stages,x_rec,y_rec,x_str,y_str,stage_plot] = stg(R_by_Rm,Rm,x_D,x_B,c_q,m_q,y,Y,eta)
    R = R_by_Rm*Rm;

    %Rectifying line
    m_r = R/(R + 1);
    c_r = x_D - m_r*x_D; %Intercept of the R line

    %Intersection points for stripping section
    x_str = (c_q - c_r)/(m_r - m_q);
    x_rec = x_str;
    y_str = m_r*x_str + c_r;
    y_rec = y_str;
    m_str = (y_str - x_B)/(x_str - x_B);
    c_str = x_B - m_str*x_B;

    %Minimum stage with efficiency added
    stages = 0;
    x_e = x_D;
    y_e = x_e;
    stage_plot = [];

    while (x_e > x_B)
        stage_plot = [stage_plot; [x_e, y_e]];
        y_temp = stage_plot(end, 2);
        x_temp = xeq(y_temp,y,Y);
        x_e = x_e - eta*(x_e - x_temp);
        stages = stages + 1;
        stage_plot = [stage_plot; [x_e, y_e]];
        y_e = min(m_r*x_e + c_r, m_str*x_e + c_str);
    end
end