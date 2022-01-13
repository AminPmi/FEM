clc
clear
close
%% -------problem 9.2 ( A first course in finite element ) ----------
% https://github.com/AminPmi
% Date: sunday, December 27, 2021
%% ------------------ Read Data From Excel File ---------------------
NodesCoordinate = xlsread('NodesCoor90');
Elements = xlsread('Elements70');
NumberOfElements = size(Elements,1);
NumberOfNodes = size(NodesCoordinate,1);

e=1:1:NumberOfElements;
node1 = NodesCoordinate(Elements(e,1),:);
node2 = NodesCoordinate(Elements(e,2),:);
node3 = NodesCoordinate(Elements(e,3),:);
node4 = NodesCoordinate(Elements(e,4),:);
CoordinateOfElement = [node1(e,:);node2(e,:);node3(e,:);node4(e,:)];

syms etha zita
nu = 0.3;  % Poisson’s ratio
E = 3e7;   % Young’s modulus
D = (E/(1-nu^2)).*[1 nu 0;nu 1 0; 0 0 (1-nu)/2];
% Shape functions
N1=(1/4)*(1-zita)*(1-etha);
N2=(1/4)*(1+zita)*(1-etha);
N3=(1/4)*(1+zita)*(1+etha);
N4=(1/4)*(1-zita)*(1+etha);
N=[N1 0 N2 0 N3 0 N4 0; 0 N1 0 N2 0 N3 0 N4];
% Gauss quadrature
zita = [ -(1/sqrt(3)); (1/sqrt(3))];
etha = [ -(1/sqrt(3)) ;(1/sqrt(3))];

% Finding K (stiffness matrix) for each element
for e=1:1:NumberOfElements
    K=zeros(8);
    for i=1:1:2
        for j=1:2
            J1=(1/4)*[etha(j)-1 1-etha(j) 1+etha(j) -etha(j)-1;zita(i)-1 -zita(i)-1 1+zita(i) 1-zita(i)]*[node1(e,:);node2(e,:);node3(e,:);node4(e,:)];
            B2=inv(J1)*(1/4)*[etha(j)-1 1-etha(j) 1+etha(j) -etha(j)-1;zita(i)-1 -zita(i)-1 1+zita(i) 1-zita(i)];
            B22=[B2(1,1) 0 B2(1,2) 0 B2(1,3) 0 B2(1,4) 0; 0 B2(2,1) 0 B2(2,2) 0 B2(2,3) 0 B2(2,4); B2(2,1) B2(1,1) B2(2,2) B2(1,2) B2(2,3) B2(1,3) B2(2,4) B2(1,4)];
            K=K+(det(J1).*(B22'*D*B22));
        end
    end
    K=mat2cell(K,8,8);
    KK(e)=K;
end
% Direct Assemblying K
K_assemble = zeros(2*NumberOfNodes);
for e=1:1:NumberOfElements
    KA{e} = zeros(2*NumberOfNodes);
    KA{e}([2*Elements(e,1)-1 2*Elements(e,1) 2*Elements(e,2)-1 2*Elements(e,2) 2*Elements(e,3)-1 2*Elements(e,3) 2*Elements(e,4)-1 2*Elements(e,4)],[2*Elements(e,1)-1 2*Elements(e,1) 2*Elements(e,2)-1 2*Elements(e,2) 2*Elements(e,3)-1 2*Elements(e,3) 2*Elements(e,4)-1 2*Elements(e,4)])=cell2mat(KK(e));
    KA;
end
K_assemble = zeros(2*NumberOfNodes);
for e=1:1:NumberOfElements
    K_assemble = K_assemble + cell2mat(KA(e));
end
% Boundary forces
syms etha zita
for e=1:1:NumberOfElements
    F_boundry = zeros(8,1);
    if e>=57
        t=[0;-(40/15)];
    else
        t=[0;0];
    end
    NN=subs(N,etha,1);
    F_boundry = F_boundry + int(NN',zita,-1,1)*t;

    F_boundry = mat2cell(double(F_boundry),8,1);
    FF(e) = F_boundry;
end
FA = zeros(2*NumberOfNodes,1);
for e=1:1:NumberOfElements
FA([2*Elements(e,1)-1 2*Elements(e,1) 2*Elements(e,2)-1 2*Elements(e,2) 2*Elements(e,3)-1 2*Elements(e,3) 2*Elements(e,4)-1 2*Elements(e,4)],1)=cell2mat(FF(e));
FA;
end
FA
% essensial DOFs
prescribedDOF=[1 2 31 32 61 62 91 92 121 122 151 152];
K_assemble_1=K_assemble;
FA1=FA;
K_assemble_1([1 2 31 32 61 62 91 92 121 122 151 152],:)=[];
K_assemble_1(:,[1 2 31 32 61 62 91 92 121 122 151 152])=[];
FA1([1 2 31 32 61 62 91 92 121 122 151 152],:)=[];
U=zeros(2*NumberOfNodes,1);
activeDof=setdiff([1:2*NumberOfNodes], [prescribedDOF]);
U1=K_assemble_1\FA1;
U(activeDof)=U1;          % Displacements

zita = [ -(1/sqrt(3)); (1/sqrt(3))];
etha = [ -(1/sqrt(3)) ;(1/sqrt(3))];

% Strains and Tensions
Strain_total=cell(NumberOfElements,1);
Sigma_total=cell(NumberOfElements,1);
for e=1:1:NumberOfElements
    for i=1:1:2
        for j=1:2
            Strain = zeros(3,1);
            Sigma = zeros(3,1);
            J11=(1/4)*[etha(j)-1 1-etha(j) 1+etha(j) -etha(j)-1;zita(i)-1 -zita(i)-1 1+zita(i) 1-zita(i)]*[node1(e,:);node2(e,:);node3(e,:);node4(e,:)];
            B22=inv(J11)*(1/4)*[etha(j)-1 1-etha(j) 1+etha(j) -etha(j)-1;zita(i)-1 -zita(i)-1 1+zita(i) 1-zita(i)];
            B222=[B22(1,1) 0 B22(1,2) 0 B22(1,3) 0 B22(1,4) 0; 0 B22(2,1) 0 B22(2,2) 0 B22(2,3) 0 B22(2,4); B22(2,1) B22(1,1) B22(2,2) B22(1,2) B22(2,3) B22(1,3) B22(2,4) B22(1,4)];
            B222=mat2cell(B222,3,8);
        end
    end
            Strain = Strain + cell2mat(B222)*U([2*Elements(e,1)-1 2*Elements(e,1) 2*Elements(e,2)-1 2*Elements(e,2) 2*Elements(e,3)-1 2*Elements(e,3) 2*Elements(e,4)-1 2*Elements(e,4)],1);
            Sigma =  Sigma + D*Strain;
            
            Strain_total{e} = Strain'; % transpose each cell
            Sigma_total{e} = Sigma';   % transpose each cell

end

%% -------------------------Plots-------------------------- %%
UU = mat2cell(U ,2*ones(size(NodesCoordinate,1),1))  % displacement of each node
for i = 1:size(NodesCoordinate,1)
    Nodes_deformed(i,1) = NodesCoordinate(i,1)+10^4*UU{i}(1);
    Nodes_deformed(i,2) = NodesCoordinate(i,2)+10^4*UU{i}(2);
end

figure('Name','Initial & Deformed Configuration','NumberTitle','off');

trisurf(Elements,NodesCoordinate(:,1),NodesCoordinate(:,2),zeros(size(NodesCoordinate,1),1),'edgecolor','b','FaceAlpha',0.1)
view(2);
set(gcf,'Color','w')
xlim([0 3])
ylim([-0.5 2])
axis equal

hold on
trisurf(Elements,Nodes_deformed(:,1),Nodes_deformed(:,2),zeros(size(NodesCoordinate,1),1),'edgecolor','b','FaceAlpha',0.4)

view(2);
set(gcf,'Color','w')
xlim([0 3])
ylim([-0.5 2])
axis equal
hold on
X=[-0.1 -0.1 0 0];
Y=[0 1 1 0];
fill(X,Y,'k');

Strain_X=zeros(size(NodesCoordinate,1),1);
Strain_total_X = cell2mat(Strain_total);
for i = 1:size(NodesCoordinate,1)
    
    Strain_X(i) = mean(Strain_total_X([find(Elements(:,1)==i);find(Elements(:,2)==i);find(Elements(:,3)==i);find(Elements(:,4)==i)],1));
    
end

figure('Name','Strain in x','NumberTitle','off');
trisurf(Elements,NodesCoordinate(:,1),NodesCoordinate(:,2),Strain_X,'edgecolor','r','facecolor','interp')
colormap jet;
set(gcf,'Color','w');
c=colorbar;
title('Strain in x');
view(2);
xlim([0 2.5])
ylim([-0.5 2.5])

Strain_Y=zeros(size(NodesCoordinate,1),1);
Strain_total_Y = cell2mat(Strain_total);
for i = 1:size(NodesCoordinate,1)
    
    Strain_Y(i) = mean(Strain_total_Y([find(Elements(:,1)==i);find(Elements(:,2)==i);find(Elements(:,3)==i);find(Elements(:,4)==i)],2));
    
end

figure('Name','Strain in y','NumberTitle','off');
trisurf(Elements,NodesCoordinate(:,1),NodesCoordinate(:,2),Strain_Y,'edgecolor','r','facecolor','interp')
colormap jet;
set(gcf,'Color','w');
c=colorbar;
title('Strain in y');
view(2);
xlim([0 2.5])
ylim([-0.5 2.5])

Stress_X=zeros(size(NodesCoordinate,1),1);
Sigma_total_X = cell2mat(Sigma_total);
for i=1:size(NodesCoordinate,1)
    
    Stress_X(i) = mean(Sigma_total_X([find(Elements(:,1)==i);find(Elements(:,2)==i);find(Elements(:,3)==i);find(Elements(:,4)==i)],1));
    
end

figure('Name','Stress in x','NumberTitle','off');
trisurf(Elements,NodesCoordinate(:,1),NodesCoordinate(:,2),Stress_X,'edgecolor','r','facecolor','interp')
colormap jet;
set(gcf,'Color','w');
c=colorbar;
title('Stress in x');
view(2);
xlim([0 2.5])
ylim([-0.5 2.5])

Stress_Y=zeros(size(NodesCoordinate,1),1);
Sigma_total_Y = cell2mat(Sigma_total);
for i=1:size(NodesCoordinate,1)
    
    Stress_Y(i) = mean(Sigma_total_Y([find(Elements(:,1)==i);find(Elements(:,2)==i);find(Elements(:,3)==i);find(Elements(:,4)==i)],2));
    
end

figure('Name','Stress in y','NumberTitle','off');
trisurf(Elements,NodesCoordinate(:,1),NodesCoordinate(:,2),Stress_Y,'edgecolor','r','facecolor','interp')
colormap jet;
set(gcf,'Color','w');
c=colorbar;
title('Stress in y');
view(2);
xlim([0 2.5])
ylim([-0.5 2.5])