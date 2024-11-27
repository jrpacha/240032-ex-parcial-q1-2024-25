clearvars 
close all
clc

numPreg = 2;

format short e
format compact

file = ['sol-prob1-q',num2str(numPreg),'.txt']; fp = fopen(file, 'w');
fprintf(fp,...
    ['%s\n\n',...
     'MN P1 ExParcial 1Q 2024-25\n',...
     'Question %d\n\n'], file, numPreg);

%Part a.1
alphasPointP = [0, 1/4, 3/4];
alphasPointQ = [1/4, 0, 3/4];
alphasPointR = [1/2, 1/2, 0];

Tp = 10.790; Tq = 30.319; Tr = 90.534;

tempAtPoints = [Tp; Tq; Tr];

A = [alphasPointP; alphasPointQ; alphasPointR];

tempAtVertexs = A\tempAtPoints;
%Hint: sum the temperatures at the vertices
sumTempAtVertexs = sum(tempAtVertexs);
fprintf(fp, ...
    ['(a.1) Temperature at the vertices: T1 = %.4e, T2 = %.4e, T3 = %.4e\n',...
     '      Hint a.1.: The sum  of the temperatures at the vertices is T1 + T2 + T3 = %.4e\n'],...
    tempAtVertexs, sumTempAtVertexs);

%Part a.2
%Interpolated temperature at the barycenter
interpTemptAtBarycenter = sumTempAtVertexs/3;
fprintf(fp,...
    '(a.2) The interpolated temperature at the barycenter is (T1 + T2 + T3)/3 = %.4e\n\n',...
    interpTemptAtBarycenter);

%Part b.1
syms x
Psi1 = expand(-x*(x-1)*(x-2)/6);
derPsi1 = diff(Psi1, x);
coefficients = coeffs(derPsi1);
coefficients = fliplr(coefficients);
fprintf(fp,...
    ['(b.1)  Psi1 = %s\n',...
     '      DPsi1 = %s\n',...
     '          a = %s\n',...
     '      Hint b.1: the apropriate shape function is of the form a·x^2 + x - 1/3\n'],...
    Psi1, derPsi1, coefficients(1));

%Part b.2
K11 = int(derPsi1*derPsi1,x,-1,2);
k11 = double(K11);
fprintf(fp,...
    '(b.2) K(k,1;1,1) = %s %s %.4e', K11, char(8776), k11);
fclose(fp);
type(file)









