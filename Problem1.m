%CHEME 5440 Final Exam Question 1
clear all;
close all;

%Set XTot, xS, and x for all solutions
xTot = [.2; .6; 1; 1.4; 1.8; 2];
xS = [0:.002:.2; 0:.006:.6; 0:.01:1; 0:.014:1.4; 0:.018:1.8; 0:.02:2];
x = xTot - xS;

%Calculate and plot forward/reverse rate
fwdRate = x .* xS;
rvsRate = xS;

figure(1)
subplot(2, 1, 1)
hold on
plot1a = plot(transpose(xS), transpose(fwdRate), 'LineWidth', 1.2);
plot1a(7) = plot(transpose(xS(6, :)), transpose(rvsRate(6, :)), 'LineWidth', 1.2);
legend(plot1a, ".2", ".6", "1", "1.4", "1.8", "2", "Reverse")
colors1a = lines(7);
set(plot1a, {'color'}, num2cell(colors1a, 2));

subplot(2, 1, 2)
hold on
plot1b = plot(transpose(xS), transpose(fwdRate), 'LineWidth', 1.2);
plot1b(7) = plot(transpose(xS(6, :)), transpose(rvsRate(6, :)), 'LineWidth', 1.2);
legend(plot1b, ".2", ".6", "1", "1.4", "1.8", "2", "Reverse")
set(plot1b, {'color'}, num2cell(colors1a, 2));
xlim([0 .1])
ylim([0 .1])

%Calculate acculation rate for xS and x, plot rates
dxdt = xS - x.*xS;
dxsdt = x.*xS - xS;

figure(2)
hold on
plot2 = plot(transpose(xS), transpose(dxdt), 'LineWidth', 1.5);
plot2(7:12) = plot(transpose(xS), transpose(dxsdt), 'LineWidth', 1.5);
legend(plot2, "dX/dt XTot = .2", "dX/dt XTot = .6", "dX/dt XTot = 1", "dX/dt XTot = 1.4", "dX/dt XTot = 1.8", "dX/dt XTot = 2", "dX*/dt XTot = .2", "dX*/dt XTot = .6", "dX*/dt XTot = 1", "dX*/dt XTot = 1.4", "dX*/dt XTot = 1.8", "dX*/dt XTot = 2")
plot([0 10], [0 0], 'k-', 'LineWidth', 1);
xlim([0 2])

colors2 = [spring(6); winter(6)];
set(plot2, {'color'}, num2cell(colors2, 2));

%Solve the ODE system and plot output
x0(1) = 1.999; 
x0(2) = .00001;
tspan = [0.001:20];

%options = odeset('NonNegative',1);
[t_out,x_out] = ode45(@(t,x) ODE3a1(t, x), tspan, x0);
%semilogx(t_out, x_out)
figure(3)
plot(t_out, x_out(:, 2))

%Plot steady-state values for X*
XTotal = [0, .2, .6, 1, 1.4, 1.8, 2];
ssXStar = [0, 0, 0, 0, .4, .8, 1];

figure(4)
hold on
plot(XTotal, ssXStar)




