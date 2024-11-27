%
%fem1DScript.m
%
clearvars
close all
clc

fileName = 'sol-prob3-q5.txt';
fp = fopen(fileName,"w","a");

fprintf(fp,['MN P3 ExParcial 1Q 2024-25\n',...
         'Question 5\n\n']);

%Real constants
E = 28.0;                          %GPa (= 28.0 kN/mm^2) 
w = 25.0e-9;                       %kN/mm^3


format long e

%Geometry
l = 0.5e3;                         %mm
L = 1.5e3;                         %mm
H = 15.0e3;                        %mm height
r = 0.5e3;                         %mm width

%Loads
P = 36.0e-6;                       %GPa (= 24.0e-6 kN/mm^2) {24,27,30,33,36}e-6 kN/mm^2

xp = H/2;                          %Point where we want to interpolate the solution

nDiv = 55;                         %Number of divisions

nodes = linspace(0, H, nDiv+1)';   %Position of the nodes
numNodes = size(nodes, 1);         %Number of nodes
h = nodes(2) - nodes(1);           %h: lenght of the nodes

elem = [1:numNodes-1; 2:numNodes]';%Connectivity matrix  
numElem = size(elem, 1);           %Number of elements 

%Compute the nodes of the elements to which point xp belongs
nod1 = ceil(xp / h); nod2 = nod1+1;

%A(x) = alpha x + beta
%A(0) = beta = rL
%A(H) = alpha*H + rL = rl <==> alpha = r*(L-l)/H

alpha = r*(l-L)/H;
beta = r*L;

stiffnessMatrix = @(x1, x2) (beta + 0.5*alpha*(x1+x2))*[1 -1; -1 1]/(x2-x1);
forceVector = @(x1,x2) (x2-x1)*(beta*[1;1]/2 + alpha*[2*x1+x2; x1+2*x2]/6);

%(a) Càlcul de les u's
K = zeros(numNodes);
F = zeros(numNodes,1);
Q = zeros(numNodes,1);
u = zeros(numNodes,1);

for e = 1:numElem

    rows = [elem(e,1); elem(e,2)];
    cols = rows;

    x1 = nodes(rows(1)); x2 = nodes(rows(2));

    Ke = E*stiffnessMatrix(x1,x2);
    K(rows,cols) = K(rows,cols) + Ke;

    Fe = forceVector(x1,x2);
    F(rows) = F(rows) + Fe;

end

%Boundary conditions B.C.
fixedNodes = 1;
freeNodes = setdiff(1:numNodes, fixedNodes);
%Natural B.C.
Q(numNodes) = -P*r*l; 

%Essential B.C.
u(1) = 0.0;

%Reduced system
Fm =  -w*F(freeNodes) + ...
    Q(freeNodes) - K(freeNodes, fixedNodes) * u(fixedNodes);
Km =  K(freeNodes, freeNodes);

%Solution of the reduced system
um = Km\Fm;
u(freeNodes) = um; solArray = [(1:numNodes)', nodes, u];

format short e
format compact

%table with the displacement at the nodes
solTable = array2table(solArray, 'VariableNames',{'Nod.','X','U'});

str = formattedDisplayText(solTable);
str = erase(str,"<strong>");
str = erase(str,"</strong>");
fprintf(fp,"%s",str);

%fprintf('%8s%13s\n','X','U')
%fprintf('%12.4e %12.4e\n',[nodes, u]')

%Plot displacements
plot(nodes,u,'o-')
xlabel('$X$ (mm)','Interpreter','latex','FontSize',15)
ylabel('$U$ (mm)','Interpreter','latex','FontSize',15)

%Part (a) Avergaged value of u
avgdU = sum(u)/numNodes;
fprintf(fp,'\nSolutions\n');
fprintf(fp,...
    '(a) <u> %s %10.4e mm. Hint. u(%.2f) %s U(%d) = %10.4e mm\n',...
    char(8776), avgdU, H, char(8776), numNodes, u(end));

%(b) Interpolated value at the midpoint of the pillar
interpU = interp1(nodes, u, xp);
fprintf(fp,...
    '(b) Interpolated value of u at x = %.2f (mm): u(%.2f) %s %10.4e mm\n',...
    xp, xp, char(8776), interpU);

%(c) Value of w s.t. the interpolated value of u at xp 3/4 the 
% value of the displacement at the top end node.
clear um;

numNodesPlusOne = numNodes + 1;
alpha = (xp-nodes(nod1))/(nodes(nod2)- nodes(nod1));

S = zeros(numNodesPlusOne);    %Extended stiffness matrix
R = zeros(numNodesPlusOne,1);  %Extended load vector
v = zeros(numNodesPlusOne,1);  %Extended solution vector, last component
                               %is the required specific weight

S(1:numNodes,1:numNodes) = K;
S(1:numNodes,numNodesPlusOne) = F;

S(numNodesPlusOne,nod1) = 1-alpha;
S(numNodesPlusOne,nod2) = alpha; 
S(numNodesPlusOne,numNodes) = -1/2;

fixedNodes = 1;
freeNodes = setdiff(1:numNodesPlusOne, fixedNodes);

%Boundary Conditions (B.C.)

%Natural B.C.
R(1:numNodes) = Q;

%Essential B.C.
v(1) = 0.0;

%Reduced system
Rm = R(freeNodes) - S(freeNodes, fixedNodes)*v(fixedNodes);
Sm = S(freeNodes, freeNodes);

%Solve the reduced system
um = Sm\Rm;
v(freeNodes) = um;

specificWeight = v(end); %Value of the required specific weight

%Check the obtained result
clear v S Sm R Rm um numNodesPlusOne alpha;

uNew = zeros(numNodes, 1); %Allocate the vector for the new displacements,
                           %i.e., the displacements corresponding to the 
                           %new spacific weight.

fixedNodes = 1;
freeNodes = setdiff(1:numNodes, fixedNodes);

%Set the essential B.C. to the new vector of the displacements, uNew.
%Note that neither the stiffness matrix K, nor the vector Q change.
uNew(1) = 0.0;

Fm =  -specificWeight*F(freeNodes) + ...
    Q(freeNodes) - K(freeNodes, fixedNodes) * uNew(fixedNodes);

%Solution of the reduced system (with the new specific weight)
um = Km\Fm;

uNew(freeNodes) = um;

% Compute the error as (the absolute value of) the difference 
% between the interpolated value of u at the midpoint and the
% 50% cent of u(numNodes)
errSpecificWeight = abs(interp1(nodes, uNew, xp) - 0.5*uNew(numNodes));

fprintf(fp,['(c) Required w = %13.6e kN/mm^3,\n\t', ...
         '   err := ||U(midpoint) - 0.5 U(top)|| = %8.2e mm\n'], ...
        specificWeight,errSpecificWeight);

fclose(fp);

type(fileName);