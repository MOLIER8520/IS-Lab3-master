clc;
clear all;
close all;

%Sukuriami duomenys 
x= 0.1:1/22:1;
y = ((1+0.6*sin(2*pi*x/0.7))+0.3*sin(2*pi*x))/2;
%I koki iejima is kokio isejimo 

% Plot original data
figure(1);
plot(x, y, 'b', 'LineWidth', 1.5);
grid on;
hold on;
title('SBF tinklo aproksimavimas');

% Aprašome duomenys SBF tinklui
c1 = 0.185; % 1 apskritimo vidurio taškas
c2 = 0.915; % 2 apskritimo vidurio taškas
r1 = 0.22; % 1 apskritimo spindulys
r2 = 0.25; % 2 apskritimo spindulys

% Sugeneruojame svorius
w1 = rand(1); 
w2 = rand(1); 
b = rand(1);
eta = 0.15; % Klaidos žingsnis

for epoch = 1:50000
    for i = 1:length(x)
        
        %skaičiuojame funkciju išėjimus
        phi1 = exp(-(x(i) - c1)^2 / (2 * r1^2));
        phi2 = exp(-(x(i) - c2)^2 / (2 * r2^2));
        
        % tinklo išėjimas
        y_hat = w1 * phi1 + w2 * phi2 + b;
        
        % klaidos skaičiavimas
        e = y(i) - y_hat;
        
        % Atnaujiname svorius
        w1 = w1 + eta * e * phi1;
        w2 = w2 + eta * e * phi2;
        b = b + eta * e;
    end
end


x_test = 0.1:1/220:1; 
Y_test = zeros(1, length(x_test));

for i = 1:length(x_test)
    
    phi1 = exp(-(x_test(i) - c1)^2 / (2 * r1^2));
    phi2 = exp(-(x_test(i) - c2)^2 / (2 * r2^2));
   
    Y_test(i) = w1 * phi1 + w2 * phi2 + b;
end

% Plot the RBF network approximation
plot(x_test, Y_test, 'r--', 'LineWidth', 1.5);
legend('Pradinis grafikas', 'SBF tinklo aproksimavimas');
hold off;
