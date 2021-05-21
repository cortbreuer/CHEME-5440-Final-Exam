%CHEME 5440 Final Exam Problem 4
clear all;
close all;

%Solve and plot solutions for the system with varying kd

%Initial conditions
x0(1) = 10; 
x0(2) = 10;
tspan = [0:.1:20];

%Set kd and solve with ODE45
kd = 2;
x0(1) = .01; 
x0(2) = .01;
[t_out,x_out1] = ode45(@(t,x) ODE4a(t, x, kd), tspan, x0);

kd = 2;
x0(1) = .1; 
x0(2) = -.1;
[t_out,x_out2] = ode45(@(t,x) ODE4a(t, x, kd), tspan, x0);

kd = 2;
x0(1) = -2; 
x0(2) = 1;
[t_out,x_out3] = ode45(@(t,x) ODE4a(t, x, kd), tspan, x0);

kd = 2;
x0(1) = 10; 
x0(2) = 10;
[t_out,x_out4] = ode45(@(t,x) ODE4a(t, x, kd), tspan, x0);

%kd = 0;
%[t_out,x_out5] = ode45(@(t,x) ODE4a(t, x, kd), tspan, x0);

%Plot all solutions
close all
figure(1)
hold on
plot(t_out, x_out1(:, 1))
plot(t_out, x_out2(:, 1))
plot(t_out, x_out3(:, 1))
plot(t_out, x_out4(:, 1))
%plot(t_out, x_out5)
legend("CA, .01", "CA, .1", "CA, 1", "CA, 10")

figure(2)
hold on
plot(t_out, x_out1(:, 2))
plot(t_out, x_out2(:, 2))
plot(t_out, x_out3(:, 2))
plot(t_out, x_out4(:, 2))
%plot(t_out, x_out5)
legend("CB, .01", "CB, .1", "CB, 1", "CB, 10")

%Phase portrait analysis

%Calculate phase vectors
size_test = 20;
R1_test = linspace(0,3,size_test);
R2_test = linspace(0,3,size_test);
kd = 4;

for i = 1:size_test
    for j = 1:size_test
        dR1(j,i) = 2 - (R1_test(i) .* R2_test(j)) - (kd .* R1_test(i));
        dR2(j,i) = (R1_test(i) .* R2_test(j)) - R2_test(j);
    end
end

%Generate quiver plot
dR1 = dR1./max(max(dR1));
dR2 = dR2./max(max(dR2));
figure(2);
plot4 = quiver(R1_test,R2_test,dR1,dR2,1.5);
axis([0 3 0 3])

%Compute the magnitude of the vectors
mags = sqrt(sum(cat(2, plot4.UData(:), plot4.VData(:), ...
            reshape(plot4.WData, numel(plot4.UData), [])).^2, 2));

%Get the current colormap
currentColormap = colormap(gca);

%Determine the color to make each arrow using a colormap
[~, ~, ind] = histcounts(mags, size(currentColormap, 1));

%Map this to a colormap to get RGB
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

%Repeat each color 3 times (using 1:3 below)
set(plot4.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');

%Repeat each color 2 times (using 1:2 below)
set(plot4.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');

%Jacobian analysis of steady-states
J = zeros(2,2);

%Jacobian for kd = 0 with one steady-state
sStateCA = 1;
sStateCB = 2;
kd = 0;
gamma = 1;

J(1, 1) =  -sStateCB - kd;
J(2, 1) = sStateCB;
J(1, 2) = sStateCA;
J(2, 2) = sStateCA - gamma;

eig0 = eig(J);
det0 = det(J);
TrJ0 = J(1,1) + J(2,2);

%Jacobian for kd = 1 with 1, 1 steady-state
sStateCA = 1;
sStateCB = 1;
kd = 1;
gamma = 1;

J(1, 1) =  -sStateCB - kd;
J(2, 1) = sStateCB;
J(1, 2) = sStateCA;
J(2, 2) = sStateCA - gamma;

eig1stab = eig(J);
det1stab = det(J);
TrJ1stab = J(1,1) + J(2,2);

%Jacobian for kd = 1 with 1, 0 steady-state
sStateCA = 1;
sStateCB = 0;
kd = 1;
gamma = 1;

J(1, 1) =  -sStateCB - kd;
J(2, 1) = sStateCB;
J(1, 2) = sStateCA;
J(2, 2) = sStateCA - gamma;

eig1un = eig(J);
det1un = det(J);
TrJ1un = J(1,1) + J(2,2);








