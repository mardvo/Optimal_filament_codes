    close all 
clear all

load('opt_reg1.mat')
figure

plot(x,A(end,:)','b','Linewidth',4)
hold on 
%plot(x,A(1,:)','b--','Linewidth',2)

F_prop
xlabel('$x$','Interpreter','Latex')
ylabel('$A$','Interpreter','Latex')
set(gca,'Fontsize',16)

t=linspace(0,2*pi/omega,15);

for i=1:15
    Y(:,i)=y2(1,:)*cos(omega*t(i))+y2(5,:)*sin(omega*t(i));
end
% Create a grayscale colormap
colormap(gray(256));

% Generate a vector of grayscale values for line colors
gray_colors = linspace(0.1, 0.9, 14);

figure
for i=2:14
    Y(:,i)=y2(1,:)*cos(omega*t(i))+y2(5,:)*sin(omega*t(i));
    line_color = [gray_colors(i) gray_colors(i) gray_colors(i)];
    plot(x, Y(:, i), 'Color', line_color, 'LineWidth', 2);
    hold on;

end
Y(:,1)=y2(1,:)*cos(omega*t(1))+y2(5,:)*sin(omega*t(1));
set(gca,'Fontsize',30)
xlabel('$x$','Interpreter','Latex','Fontsize',30)
ylabel('$y$','Interpreter','Latex','Fontsize',30)
plot(x, Y(:, 1), 'k-','LineWidth', 3);
% Adjust the colormap range to grayscale
colormap gray;

% Hold off to prevent further plots from overwriting the current one
hold off;
 