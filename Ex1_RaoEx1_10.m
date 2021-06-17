%Solves Example 1.10 of (Rao - old edition)
%This code is edited By Fizza Shamim & Ali Irshayyid
% The problem description is as follows.
% Minimize f1(X)=p*l*x1*x2 and
% f2(X)=-sqrt(E*x1*(x2^3))/(4*(l^2)*(M+(33/140)*f1))
% subject to:
% g1(X): (M*9.8)/(x(1)*x(2))-sigma<=0
% g2(X): (M*9.8)/(x(1)*x(2))-(((pi)^2*E*x(2)^2)/(48*l^2))<=0
% g3(X): -x(1)<=0
% g4(x): -x(2)<=0
% g5(x): -W <= 10
% g6(x): W <= 100

clc; clear all;
state = 0;
input_type = 0;
global p E M max_sigma lambda l
M=37854.117;
%% Asking the user for information
while (state ==0)
fprintf('HELLO USER!\nWELCOME TO EX1.10\nKINDLY, ENTER THE MATERIAL TYPE \n');
s1=input('concrete or steel:       ','s');
while(input_type == 0)
    
l=input('KINDLY, ENTER THE COLUMN LENGTH YOU WOULD PREFFER IN METERS:    ');
if ( isnumeric(str2double(l)))
    input_type =1;
else
    fprintf("Invalid length input, Please try again \n")
end
end
%% Assigning Values based on the user input
if (strcmp('concrete',s1)==1)
    p=2398; %Density, kg/m^3
    E=30;   %Young's Modulus, GPA
    max_sigma=40; %Maximum Compressive stress
    state = 1;
elseif (strcmp('steel',s1)==1)
  p=0.00785; %Density, kg/m^3
    E=200;   %Young's Modulus, GPA
    max_sigma=400; %Maximum Compressive stress
    state = 1;
else 
    fprintf("Invalid Input material, Please try again \n \n \n")
    
    
end
end

%% Implementing the optimization
lambda_array = (0.01:0.05:1);
x_min = zeros(3,20);
f_min = zeros(1,20);
for i = 1:1:20

lambda = lambda_array(i);
%Calls the functions: ex1_4fun and ex1_4con
x0=[7; 1000];
lb=[0; 0];
[x, funcion_value]=fmincon(@mini_fun,x0,[],[],[],[],lb,[],@mini_constraints);
f_min(1,i) = funcion_value;
g = mini_constraints(x);
x_min(1,i) = x(1);
x_min(2,i) = x(2);

if (g(1) <= 0 && g(2) <=0)
    disp("this min is valid");
    x_min(3,i) = 1;
else 
    disp("This X_min doesn't satisfy the constraints !!!!");
end
end
%% Plotting 
plot (lambda_array,f_min );
title("Fmin changing w.r.t lambda");
xlabel("Lambda");
ylabel("Fmin");

%% Functions
function f=mini_fun(X)
%f=0.7402*(x(1)^2)*x(2)*x(3);
% M=37534.6; l=2;  p=2398; %Density, kg/m^3
% E=30; max_sigma=40;
global p E M l lambda
x1=X(1);
x2=X(2);

f1=p*l*x1*x2;
f_2=(E*x1*(x2^3))/(4*(l^2)*(M+(33/140)*f1));
f2=-sqrt(f_2);
f= lambda * f1+(1 - lambda)* f2;
% f= f1+f2;
end 

function [g,geq]=mini_constraints(X)
% M=37534.6; l=2;  p=2398; %Density, kg/m^3
% E=30; max_sigma=40;
% load('workspace');
global  E M max_sigma l p
x1=X(1);
x2=X(2);
z=(M*9.8)/(x1*x2);
g(1)=z-max_sigma;
g(2)=z-(((pi^2)*E*x2^2)/(48*l^2));

% strict the Frequency of Transversal Vibration approximately between 10 to 100 rad/sec.
w=(E*x1*(x2^3))/(4*(l^2)*(M+(33/140)*p*l*x1*x2));
g(3)=10 - sqrt(w); % W >= 10
g(4)= sqrt(w) - 100; % W <= 100
geq=[];
end

