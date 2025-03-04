clear all;
close all;

npts=1000;

x = linspace(0,1,npts);
A0=0.07931; % maximum value of A

consts; % solve main problem
const_cong; % solve the conjugate problem

%dF/dA
figure
plot(x,-dFdA,'LineWidth',2)
xlabel('$x$','Interpreter','Latex')
ylabel('$dF/dA$','Interpreter','Latex')
set(gca,'Fontsize',16)


% figure for time snaps
t=linspace(0,2*pi,15);

% Create a grayscale colormap
colormap(gray(256));

% Generate a vector of grayscale values for line colors
gray_colors = linspace(0.1, 0.9, 14);

figure
for i=2:14
    Y(:,i)=(real(y)*cos(t(i))-imag(y)*sin(t(i)));
    line_color = [gray_colors(i) gray_colors(i) gray_colors(i)];
    plot(x, Y(:, i), 'Color', line_color, 'LineWidth', 2);
    hold on;

end

Y(:,1)=real(y)*cos(t(1))-imag(y)*sin(t(1));
set(gca,'Fontsize',30)
xlabel('$x$','Interpreter','Latex','Fontsize',30)
ylabel('$y$','Interpreter','Latex','Fontsize',30)
plot(x, Y(:, 1), 'k-','LineWidth', 3);
% Adjust the colormap range to grayscale
colormap gray;
hold off;

