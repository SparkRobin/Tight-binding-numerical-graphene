% 双层石墨烯
% Generate twist bilayer graphene LAMMPS file 
% Output is atoms in a supercell 用ovito打开

% Energy spectrum and quantum Hall effect in twisted bilayer graphene
% Pilkyung Moon and Mikito Koshino
% Phys. Rev. B 85, 195458 – Published 29 May 2012

clear
close all
Grab=zeros(100,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=22
n=m+1
theta=acos(0.5*(m^2+n^2+4*m*n)/(m^2+n^2+m*n)); %扭转角度

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_deg=theta/pi*180
acc = 1.42; % C-C bond length
lattice_constant = acc * sqrt(3); % lattice constant

dx = lattice_constant;
dy = lattice_constant*sqrt(3)/2;

ni = 3*max([m,n]); 
nj = ni;
Gra=[];


ii= 0;
for i=-10:ni
    for j=-10:nj
        ii=ii+1;
        if mod(j,2)==1
            Gra(ii,1)=(i-1)*dx;
            Gra(ii,2)=(j-1)*dy;
        else
            Gra(ii,1)=(i-1)*dx+dx/2;  %construct triangle lattice
            Gra(ii,2)=(j-1)*dy;
        end
    end
end

nB=ii;

for i=1:nB
    Gra(i+nB,1)=Gra(i,1);
    Gra(i+nB,2)=Gra(i,2)+acc;
end

% Gra=(RotateZ(pi/6)*Gra')';


Grat=(RotateZ(theta)*Gra')';

%Calculate supercell boundary
a1=[lattice_constant 0];
a2=[1/2*lattice_constant sqrt(3)*lattice_constant/2];
L=zeros(4,2);
% Cut boundary
L(2,:)=[m*a1(1)+n*a2(1) n*a2(2)];
L(4,:)=[L(2,1)*cos(pi/3)-L(2,2)*sin(pi/3) L(2,1)*sin(pi/3)+L(2,2)*cos(pi/3)];
L(3,:)=[L(2,1)+L(4,1) L(2,2)+L(4,2)];

[GraBin,on]=inpolygon(Gra(:,1),Gra(:,2),L(:,1),L(:,2)); 
GraBin=GraBin|on;
[GraTin,on]=inpolygon(Grat(:,1),Grat(:,2),L(:,1),L(:,2)); 
GraTin=GraTin|on;
GraB=Gra(GraBin,:);GraT=Grat(GraTin,:);


% Make x axis parallel to L1

kk=-atan(L(2,2)/L(2,1));  %angle between L1 and x axis

GraB=(RotateZ(kk)*GraB')';
GraT=(RotateZ(kk)*GraT')';
GraB(:,3)=0;
GraT(:,3)=3.35;

L=(RotateZ(kk)*L')';

Length=L(2,1);
% %delete some atoms in the corner to ensure not duplicated 
Gra=[GraB;GraT];
%
xl=Length;
yl=Length*sqrt(3)/2;
% Gra(Gra(:,1)==max(Gra(:,1)),:)=[];
nnk1=[];nnk2=[];nnk3=[];
for i=1:length(Gra)
    if abs(Gra(i,1)-xl)<1e-2 && Gra(i,2)<1e-2
        nnk1=[nnk1 i];
    end
    if abs(Gra(i,1)-3*xl/2)<1e-2
        nnk2=[nnk2 i];
    end
    if abs(Gra(i,1)-xl/2)<1e-2 && abs(Gra(i,2)-yl)<1e-2
        nnk3=[nnk3 i];
    end
end
Gra([nnk1 nnk2 nnk3],:)=[];

figure
hold on 
scatter(GraB(:,1),GraB(:,2),'blue','filled')
axis equal
scatter(GraT(:,1),GraT(:,2),'red')
plot([L(1,1),L(2,1)],[L(1,2),L(2,2)],LineWidth=2,Color='y');
plot([L(2,1),L(3,1)],[L(2,2),L(3,2)],LineWidth=2,Color='y');
plot([L(3,1),L(4,1)],[L(3,2),L(4,2)],LineWidth=2,Color='y');
plot([L(4,1),L(1,1)],[L(4,2),L(1,2)],LineWidth=2,Color='y');
%% output
GraB=Gra(Gra(:,3)==0,:);
GraT=Gra(Gra(:,3)~=0,:);
nB=length(GraB);
nT=length(GraT);
iG=length(Gra)

 % boundary of the graphene
XGra=[min(Gra(:,1)),xl];
YGra=[min(Gra(:,2)),yl];
ZGra=[min(Gra(:,3)),max(Gra(:,3))];

FileName=sprintf('m%dn%dTBG.in',m,n);
fid=fopen(FileName,'wt'); %Creat a file to write
fprintf(fid,'#Lammps data file\n\n');
fprintf(fid,'%5d atoms\n',iG); % totoal atoms number
fprintf(fid,'%5d atom types\n\n',2); % types number

% boundaries
fprintf(fid,'%8.4f %8.4f xlo xhi\n',XGra(1),XGra(2));
fprintf(fid,'%8.4f %8.4f ylo yhi\n',YGra(1),YGra(2));
fprintf(fid,'%8.4f %8.4f zlo zhi\n',ZGra(1)-10,ZGra(2)+10);
fprintf(fid,'%8.4f %8.4f %8.4f xy xz yz \n\n',1/2*xl,0,0); % Tricilinic box
%mass
fprintf(fid,'Masses\n\n');
fprintf(fid,' 1 12\n'); 
fprintf(fid,' 2 12\n\n');

fprintf(fid,'Atoms\n\n');

for i=1:nT
    fprintf(fid,'%5d %1d %7.3f  %7.3f  %7.3f\n',i,1,GraT(i,:));
    
end
for i=1:nB
    fprintf(fid,'%5d %1d %7.3f  %7.3f  %7.3f\n',i+nT,2,GraB(i,:));
    
end

fclose(fid);

