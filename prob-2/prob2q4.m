clearvars
close all

numPreg = 4;

p = [16,5];

load meshMontseny.mat
numNodes = size(nodes,1);
numElem = size(elem,1);

file = ['sol-prob2-q',num2str(numPreg),'.txt']; fp = fopen(file, 'w');

fprintf(fp,...
    ['MN-P2-ExParcial-1Q-2024-25\n',...
     'Question %d\n\n',...
     'Consider the triangular mesh meshMontseny.mat that you have available in this folder.\n',...
     'The mesh is from the Montseny Natural Park, with units (the positions of the nodes) in\n',...
     'kilometers. Consider the point P = (%.1f, %.1f).\n\n'], numPreg, p);

%Part (a) (3 points)
%The mean of the 𝑥-component of the three nodes of the element containing 𝑃 is:


for e = 1:numElem
    nods = elem(e,:);
    vertexs = nodes(nods,:);
    [aphas, isInside] = baryCoord(vertexs, p);
    if isInside > 0
        meanX = sum(vertexs(:,1))/3; meanY = sum(vertexs(:,2))/3;
        fprintf(fp,...
            ['(a) (3 points) The mean of the x-component of the element containig P is:\n\n',...
             '    Sol. <x> = %.4e\n\n',...
             '    Hint. The mean of the y-component of the element containig P is: <y> = %.4e\n\n'],...
            meanX,meanY);
        break
    end
end
 %Part (b)
 elemNum = 1201;
 q = sum(nodes(elem(elemNum,:),:))/3; distPQ = norm(p-q);
 fprintf(fp,...
     ['(b) (3 points) If the point Q is the barycenter of element %d, de distance between P\n',...
      '    and Q is:\n\n',...
      '    Sol. dist(P, Q) = %.4e\n\n'], elemNum, distPQ);
%Part (c)
r = [17,8];
Tp = 18.7; Tq = 17.6; Tr = 17.2;
tempAtVertexs = [Tp; Tq; Tr];
s = [18,6];
vertexs = [p; q; r];
[alphas, isInside] = baryCoord(vertexs, s);
interpTempAtS = alphas*tempAtVertexs;
fprintf(fp,...
    ['(c) (2 points) Let us consider the point R = (%.1f,%.1f). If the temperature at points P,\n',...
     '    Q, and R is T(P) = 18.7, T(Q) = 17.6, T(R) = 17.2, what will be the interpolated\n',... 
     '    temperature at point S = (%.1f,%.1f)?\n\n',...
     '    Sol. Interpolated temperature at S: %.4e\n\n'], r, s, interpTempAtS)
%Part (d)
Surface = 315.0;
areaMesh = 0.0;
columnOfOnes = ones(3,1);
for e = 1:numElem
    nods = elem(e,:);
    vertexs = nodes(nods,:);
    areaMesh = areaMesh + det([columnOfOnes, vertexs]);
end
relErr = abs(Surface - 0.5*areaMesh)/Surface;
fprintf(fp,...
    ['(d) (2 points) In a book, we saw that the surface of the Natural Park is S = %.1f km^2.\n',... 
     '    What is the relative error of the area defined by our mesh if we take S as the true\n',...
     '    value?\n\n',...
     '    Sol. Relative Error = |Measured Value - True Value| / |True Value| = %.4e\n'],...
     Surface,relErr);

fclose(fp);
type(file);