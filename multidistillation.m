function multidistillation(~)
clc
%This function uses several different parts and ultimately outputs a 
%Tray Number vs. Composition diagram for a multi-component mixture.
%As of right now, the convergence of the program relies heavily on the
%initial information provided. This program assumes constant relative
%volatility throughout the column.

%--------------------------------------------------------------------------

%User inputted information

z_F = [.25 .25 .25 .25]; %Feed composition for each component
x_D = [.504 .477 .0186 0.00000767]; %Desired distillate composition
x_B = [0.000233 .02 .48 .50]; %Desired bottoms composition
R_by_Rm = 4/.62; %R_by_Rm is the value of R/Rm
q = 1; %Condition of the feed entering the column
l = 1; %Liquid fraction of feed entering the column
alpha = [5.0 2.5 1.0 0.2]; %Relative volatility compared to one component
F = 100; %Molar flow of feed stream

%--------------------------------------------------------------------------

%Calculates the minimum reflux using the Underwood Method required for
%distillation

x_0 = 1.1; %Initial guess for theta, which should be between values of alpha
q_fun = @(theta)(sum((alpha.*z_F)./(alpha-theta)) - (1 - q));
options = optimoptions('fsolve', 'Display', 'off');
thetaC = fsolve(q_fun, x_0, options);
Rm = sum((alpha.*x_D)./(alpha-thetaC)) - 1;
fprintf('The minimum reflux required is %.2f \n', Rm);

R = R_by_Rm * Rm;
fprintf('The reflux used is %.2f \n', R);

%--------------------------------------------------------------------------

%Calculates number of trays using the Fenske Eqs. and the feed location
%using the Kirkbride Eqs. required for distillation

B = F*(z_F(1) - x_D(1))/(x_B(1) - x_D(1));
D = F - B;
alphaL = max(alpha);
i = find(alpha==alphaL);
alphaH = min(alpha);
j = find(alpha==alphaH);
xLD = x_D(i);
xHD = x_D(j);
xHB = x_B(j);
xLB = x_B(i);
xHF = z_F(j);
xLF = z_F(i);
Ntot = log((xLD/xHD)*(xHB/xLB))/log(alphaL);

fprintf('The estimated number of stages is %.2f \n', Ntot);

x_0 = 17; %Initial guess for number of stages
n_fun = @(n)(0.206*log((xHF/xLF)*(B/D)*(xLB/xHD)^2) - log((Ntot-n)/n));
options = optimoptions('fsolve', 'Display', 'off');
Nf = fsolve(n_fun, x_0, options);

fprintf('The the feed tray is tray %.0f \n', Nf);

%--------------------------------------------------------------------------

%Determines material balances using the reflux ratio and completes the
%design using Tray-to-Tray Method

l_str = (R*D + F*l)/B;
v_str = D*(1+R)/B;
l_rec = R;
v_rec = 1+R;

x(1,:) = x_B;
s = size(z_F);
Ntot = 7;
Nf = 5;

for i = 1:Ntot
    for j = 1:s(2)
        if i <= Nf
            if i == 1
                y(i,j) = alpha(j)*x(i,j)/(sum(alpha.*x(i,:)));
            else
                x(i,j) = (v_str*y(i-1,j) + x(1,j))/l_str;
            end
        else
            x(i,j) = (v_rec*y(i-1,j) - x_D(j))/l_rec;
        end 
    end
    y(i,:) = alpha.*x(i,:)/(sum(alpha.*x(i,:))); 
end

%--------------------------------------------------------------------------

%Plotting

c = size(x);
x([c(1)],:) = y([c(1)],:);
stage_arr = [1:1:Ntot];
figure('Name','Multicomponent Distillation')
plot(x,stage_arr)
xlabel('Composistion, mole fraction');
ylabel('Tray Number');
grid on
title('Tray Composition vs. Position in Tower')

end