%Nlines = 8;
%thresold = 0.2;
%K = 10;
%VV = [9     9     5     1; 6     3     7     5; 2    10    10     3; 4     1     3     9; 9     5     4     6; 2     5     3    10; 4     5     3    10; 10     2     2     6];

close all
clear
%Nlines = 30;
%thresold = 3;
%K = 50;
%VV = randi(K,Nlines,4);

% 2D joint generator
%clear
%hold on
fillings=[];
tic
%random numbers
start1=rem(now,1);
start2=start1*1000000;
%rand('state',start2)
rng(start2,'v5uniform');
% boundary
xlo=0;
xup=30;
ylo=0;
yup=30;


JOINTS.beginnings=[];
JOINTS.endings=[];
% Number of joint sets

njointsets=input('Number of joint sets? (1-4) ');

disp ('Set 1:')
for i=1:njointsets
JOINTS.(genvarname(['Set', num2str(i)])).amean=input('Mean dip angle of joints? ');
JOINTS.(genvarname(['Set', num2str(i)])).astandard=input('Standard deviation of the dip angle? ');
JOINTS.(genvarname(['Set', num2str(i)])).lmean=input('Mean trace length? ');
%random=0.5+8*rand(1,1);
%X=['Tracelength=',num2str(random)];
%disp(X)
%JOINTS.(genvarname(['Set', num2str(i)])).lmean=random;

JOINTS.(genvarname(['Set', num2str(i)])).onedimfreq=input('Joint frequency? ');
%random=0.5+4*rand(1,1);
%X=['Frequency=',num2str(random)];
%disp(X)
%JOINTS.(genvarname(['Set', num2str(i)])).onedimfreq=random;
% number of joints
JOINTS.(genvarname(['Set', num2str(i)])).njoint=ceil(JOINTS.(genvarname(['Set', num2str(i)])).onedimfreq*xup*(yup/JOINTS.(genvarname(['Set', num2str(i)])).lmean));
JOINTS.(genvarname(['Set', num2str(i)])).angle=normrnd(JOINTS.(genvarname(['Set', num2str(i)])).amean,JOINTS.(genvarname(['Set', num2str(i)])).astandard,JOINTS.(genvarname(['Set', num2str(i)])).njoint,1);
% 2D joint point generation (Poisson process)
JOINTS.(genvarname(['Set', num2str(i)])).x1=xlo+(xup-xlo)*rand(JOINTS.(genvarname(['Set', num2str(i)])).njoint,1);
JOINTS.(genvarname(['Set', num2str(i)])).y1=ylo+(yup-ylo)*rand(JOINTS.(genvarname(['Set', num2str(i)])).njoint,1);

% dip angle of joint
% Normal distribution
%JOINTS.(genvarname(['Set', num2str(i)])).amean=input('Mean dip angle of joints? ');
%JOINTS.(genvarname(['Set', num2str(i)])).astandard=input('Standard deviation of the dip angle? ');
JOINTS.(genvarname(['Set', num2str(i)])).angle=normrnd(JOINTS.(genvarname(['Set', num2str(i)])).amean,JOINTS.(genvarname(['Set', num2str(i)])).astandard,JOINTS.(genvarname(['Set', num2str(i)])).njoint,1);

% Simulation of the trace length distribution

% uniform distribution
% lmean=input('Trace length mean')
% random1=rand(njoint,1)
% length=(lmean*2).*random1

% neg. exponential distribution

JOINTS.(genvarname(['Set', num2str(i)])).length=exprnd(JOINTS.(genvarname(['Set', num2str(i)])).lmean,JOINTS.(genvarname(['Set', num2str(i)])).njoint,1);

% log. normal distribution
% lmean=input('Trace length mean')
% lstandard=input('standard deviation')
% length=lognrnd(lmean,lstandard,njoint,1)

% calculation of x2,y2
JOINTS.(genvarname(['Set', num2str(i)])).lh=JOINTS.(genvarname(['Set', num2str(i)])).length./2; %replace with lmean with length
JOINTS.(genvarname(['Set', num2str(i)])).dx=JOINTS.(genvarname(['Set', num2str(i)])).lh.*cos(JOINTS.(genvarname(['Set', num2str(i)])).angle*pi/180);
JOINTS.(genvarname(['Set', num2str(i)])).dy=JOINTS.(genvarname(['Set', num2str(i)])).lh.*sin(JOINTS.(genvarname(['Set', num2str(i)])).angle*pi/180);

JOINTS.(genvarname(['Set', num2str(i)])).x2=JOINTS.(genvarname(['Set', num2str(i)])).x1+JOINTS.(genvarname(['Set', num2str(i)])).dx;
JOINTS.(genvarname(['Set', num2str(i)])).y2=JOINTS.(genvarname(['Set', num2str(i)])).y1+JOINTS.(genvarname(['Set', num2str(i)])).dy;

JOINTS.(genvarname(['Set', num2str(i)])).x3=JOINTS.(genvarname(['Set', num2str(i)])).x1-JOINTS.(genvarname(['Set', num2str(i)])).dx;
JOINTS.(genvarname(['Set', num2str(i)])).y3=JOINTS.(genvarname(['Set', num2str(i)])).y1-JOINTS.(genvarname(['Set', num2str(i)])).dy;


% cutting x2,y2
JOINTS.(genvarname(['Set', num2str(i)])).i1=find(JOINTS.(genvarname(['Set', num2str(i)])).x2>xup);
JOINTS.(genvarname(['Set', num2str(i)])).dxni=xup-JOINTS.(genvarname(['Set', num2str(i)])).x1(JOINTS.(genvarname(['Set', num2str(i)])).i1,1);
JOINTS.(genvarname(['Set', num2str(i)])).dyni=JOINTS.(genvarname(['Set', num2str(i)])).dxni.*tan(JOINTS.(genvarname(['Set', num2str(i)])).angle(JOINTS.(genvarname(['Set', num2str(i)])).i1,1)*pi/180);
JOINTS.(genvarname(['Set', num2str(i)])).dx(JOINTS.(genvarname(['Set', num2str(i)])).i1,1)=JOINTS.(genvarname(['Set', num2str(i)])).dxni;
JOINTS.(genvarname(['Set', num2str(i)])).dy(JOINTS.(genvarname(['Set', num2str(i)])).i1,1)=JOINTS.(genvarname(['Set', num2str(i)])).dyni;
JOINTS.(genvarname(['Set', num2str(i)])).x2=JOINTS.(genvarname(['Set', num2str(i)])).x1+JOINTS.(genvarname(['Set', num2str(i)])).dx;
JOINTS.(genvarname(['Set', num2str(i)])).y2=JOINTS.(genvarname(['Set', num2str(i)])).y1+JOINTS.(genvarname(['Set', num2str(i)])).dy;


JOINTS.(genvarname(['Set', num2str(i)])).k1=find(JOINTS.(genvarname(['Set', num2str(i)])).y2>yup);
JOINTS.(genvarname(['Set', num2str(i)])).dynk=yup-JOINTS.(genvarname(['Set', num2str(i)])).y1(JOINTS.(genvarname(['Set', num2str(i)])).k1,1);
JOINTS.(genvarname(['Set', num2str(i)])).dxnk=JOINTS.(genvarname(['Set', num2str(i)])).dynk./tan(JOINTS.(genvarname(['Set', num2str(i)])).angle(JOINTS.(genvarname(['Set', num2str(i)])).k1,1)*pi/180);
JOINTS.(genvarname(['Set', num2str(i)])).dx(JOINTS.(genvarname(['Set', num2str(i)])).k1,1)=JOINTS.(genvarname(['Set', num2str(i)])).dxnk;
JOINTS.(genvarname(['Set', num2str(i)])).dy(JOINTS.(genvarname(['Set', num2str(i)])).k1,1)=JOINTS.(genvarname(['Set', num2str(i)])).dynk;
JOINTS.(genvarname(['Set', num2str(i)])).y2=JOINTS.(genvarname(['Set', num2str(i)])).y1+JOINTS.(genvarname(['Set', num2str(i)])).dy;
JOINTS.(genvarname(['Set', num2str(i)])).x2=JOINTS.(genvarname(['Set', num2str(i)])).x1+JOINTS.(genvarname(['Set', num2str(i)])).dx;


% cutting x3,y3


JOINTS.(genvarname(['Set', num2str(i)])).n1=find(JOINTS.(genvarname(['Set', num2str(i)])).x3<xlo);
JOINTS.(genvarname(['Set', num2str(i)])).dxnn=JOINTS.(genvarname(['Set', num2str(i)])).x1(JOINTS.(genvarname(['Set', num2str(i)])).n1,1)-xlo;
JOINTS.(genvarname(['Set', num2str(i)])).dynn=JOINTS.(genvarname(['Set', num2str(i)])).dxnn.*tan(JOINTS.(genvarname(['Set', num2str(i)])).angle(JOINTS.(genvarname(['Set', num2str(i)])).n1,1)*pi/180);
JOINTS.(genvarname(['Set', num2str(i)])).dx(JOINTS.(genvarname(['Set', num2str(i)])).n1,1)=JOINTS.(genvarname(['Set', num2str(i)])).dxnn;
JOINTS.(genvarname(['Set', num2str(i)])).dy(JOINTS.(genvarname(['Set', num2str(i)])).n1,1)=JOINTS.(genvarname(['Set', num2str(i)])).dynn;
JOINTS.(genvarname(['Set', num2str(i)])).x3=JOINTS.(genvarname(['Set', num2str(i)])).x1-JOINTS.(genvarname(['Set', num2str(i)])).dx;
JOINTS.(genvarname(['Set', num2str(i)])).y3=JOINTS.(genvarname(['Set', num2str(i)])).y1-JOINTS.(genvarname(['Set', num2str(i)])).dy;

JOINTS.(genvarname(['Set', num2str(i)])).o1=find(JOINTS.(genvarname(['Set', num2str(i)])).y3<ylo);
JOINTS.(genvarname(['Set', num2str(i)])).dyno=JOINTS.(genvarname(['Set', num2str(i)])).y1(JOINTS.(genvarname(['Set', num2str(i)])).o1,1)-ylo;
JOINTS.(genvarname(['Set', num2str(i)])).dxno=JOINTS.(genvarname(['Set', num2str(i)])).dyno./tan(JOINTS.(genvarname(['Set', num2str(i)])).angle(JOINTS.(genvarname(['Set', num2str(i)])).o1,1)*pi/180);
JOINTS.(genvarname(['Set', num2str(i)])).dx(JOINTS.(genvarname(['Set', num2str(i)])).o1,1)=JOINTS.(genvarname(['Set', num2str(i)])).dxno;
JOINTS.(genvarname(['Set', num2str(i)])).dy(JOINTS.(genvarname(['Set', num2str(i)])).o1,1)=JOINTS.(genvarname(['Set', num2str(i)])).dyno;
JOINTS.(genvarname(['Set', num2str(i)])).y3=JOINTS.(genvarname(['Set', num2str(i)])).y1-JOINTS.(genvarname(['Set', num2str(i)])).dy;
JOINTS.(genvarname(['Set', num2str(i)])).x3=JOINTS.(genvarname(['Set', num2str(i)])).x1-JOINTS.(genvarname(['Set', num2str(i)])).dx;



% new length
JOINTS.(genvarname(['Set', num2str(i)])).lnew=sqrt(JOINTS.(genvarname(['Set', num2str(i)])).dx.*JOINTS.(genvarname(['Set', num2str(i)])).dx+JOINTS.(genvarname(['Set', num2str(i)])).dy.*JOINTS.(genvarname(['Set', num2str(i)])).dy);
% total joint length
JOINTS.(genvarname(['Set', num2str(i)])).lsum=sum(JOINTS.(genvarname(['Set', num2str(i)])).lnew.*2);

% plot of points
% plot (x1,y1,'b+')
% hold on
%plot (x2,y2,'r*')

JOINTS.(genvarname(['Set', num2str(i)])).crackstring=zeros(JOINTS.(genvarname(['Set', num2str(i)])).njoint,1);
JOINTS.(genvarname(['Set', num2str(i)])).crackstring(:,1)=999999;
JOINTS.(genvarname(['Set', num2str(i)])).cracks=[JOINTS.(genvarname(['Set', num2str(i)])).crackstring,JOINTS.(genvarname(['Set', num2str(i)])).x3,JOINTS.(genvarname(['Set', num2str(i)])).y3,JOINTS.(genvarname(['Set', num2str(i)])).x2,JOINTS.(genvarname(['Set', num2str(i)])).y2];

% plot of fracture map
JOINTS.(genvarname(['Set', num2str(i)])).dxneu=JOINTS.(genvarname(['Set', num2str(i)])).x2-JOINTS.(genvarname(['Set', num2str(i)])).x3;
JOINTS.(genvarname(['Set', num2str(i)])).dyneu=JOINTS.(genvarname(['Set', num2str(i)])).y2-JOINTS.(genvarname(['Set', num2str(i)])).y3;
%quiver(x3,y3,dxneu,dyneu,0)
%plot(x3,y3)
%axis([xlo,xup,ylo,yup]);
%axis square;

JOINTS.(genvarname(['Set', num2str(i)])).beginnings=[JOINTS.(genvarname(['Set', num2str(i)])).x2 JOINTS.(genvarname(['Set', num2str(i)])).x3]';
JOINTS.(genvarname(['Set', num2str(i)])).endings=[JOINTS.(genvarname(['Set', num2str(i)])).y2 JOINTS.(genvarname(['Set', num2str(i)])).y3]';

%plot(JOINTS.(genvarname(['Set', num2str(i)])).beginnings, JOINTS.(genvarname(['Set', num2str(i)])).endings,'color','k');


JOINTS.beginnings=horzcat(JOINTS.beginnings,JOINTS.(genvarname(['Set',num2str(i)])).beginnings);
JOINTS.endings=horzcat(JOINTS.endings,JOINTS.(genvarname(['Set',num2str(i)])).endings);

if i<njointsets
setdescr=['Set ' num2str(i+1) ':'];

disp(setdescr);
end
end
JOINTS.X1Y1X2Y2=(vertcat(JOINTS.beginnings,JOINTS.endings))';
if njointsets==1
    VV1(:,1)=JOINTS.Set1.x2;
    VV1(:,2)=JOINTS.Set1.y2;
    VV1(:,3)=JOINTS.Set1.x3;
    VV1(:,4)=JOINTS.Set1.y3;
    VV=vertcat(VV1);
end
if njointsets==2
    VV1(:,1)=JOINTS.Set1.x2;
    VV1(:,2)=JOINTS.Set1.y2;
    VV1(:,3)=JOINTS.Set1.x3;
    VV1(:,4)=JOINTS.Set1.y3;
    VV2(:,1)=JOINTS.Set2.x2;
    VV2(:,2)=JOINTS.Set2.y2;
    VV2(:,3)=JOINTS.Set2.x3;
    VV2(:,4)=JOINTS.Set2.y3;
	VV=vertcat(VV1,VV2);
end
if njointsets==3
    VV1(:,1)=JOINTS.Set1.x2;
    VV1(:,2)=JOINTS.Set1.y2;
    VV1(:,3)=JOINTS.Set1.x3;
    VV1(:,4)=JOINTS.Set1.y3;
    VV2(:,1)=JOINTS.Set2.x2;
    VV2(:,2)=JOINTS.Set2.y2;
    VV2(:,3)=JOINTS.Set2.x3;
    VV2(:,4)=JOINTS.Set2.y3;
    VV3(:,1)=JOINTS.Set3.x2;
    VV3(:,2)=JOINTS.Set3.y2;
    VV3(:,3)=JOINTS.Set3.x3;
    VV3(:,4)=JOINTS.Set3.y3;
	VV=vertcat(VV1,VV2,VV3);
end
if njointsets==4
    VV1(:,1)=JOINTS.Set1.x2;
    VV1(:,2)=JOINTS.Set1.y2;
    VV1(:,3)=JOINTS.Set1.x3;
    VV1(:,4)=JOINTS.Set1.y3;
    VV2(:,1)=JOINTS.Set2.x2;
    VV2(:,2)=JOINTS.Set2.y2;
    VV2(:,3)=JOINTS.Set2.x3;
    VV2(:,4)=JOINTS.Set2.y3;
    VV3(:,1)=JOINTS.Set3.x2;
    VV3(:,2)=JOINTS.Set3.y2;
    VV3(:,3)=JOINTS.Set3.x3;
    VV3(:,4)=JOINTS.Set3.y3;
    VV4(:,1)=JOINTS.Set4.x2;
    VV4(:,2)=JOINTS.Set4.y2;
    VV4(:,3)=JOINTS.Set4.x3;
    VV4(:,4)=JOINTS.Set4.y3;
	VV=vertcat(VV1,VV2,VV3,VV4);
end

beginnings=[VV(:,1),VV(:,3)]';
endings=[VV(:,2),VV(:,4)]';
bordersbeginnings=[xlo,xup;xlo,xup;xlo,xlo;xup,xup];
bordersendings=[ylo,ylo;yup,yup;ylo,yup;ylo,yup];
borders=[bordersbeginnings(:,1),bordersendings(:,1),bordersbeginnings(:,2),bordersendings(:,2)];
HORIZONTALS=[xlo,2,xup,2;xlo,4,xup,4;xlo,6,xup,6;xlo,8,xup,8;xlo,10,xup,10;xlo,12,xup,12;xlo,14,xup,14;xlo,16,xup,16;xlo,18,xup,18;xlo,20,xup,20;xlo,22,xup,22;xlo,24,xup,24;xlo,26,xup,26;xlo,28,xup,28;xlo,30,xup,30;xlo,32,xup,32;xlo,34,xup,34;xlo,36,xup,36;xlo,38,xup,38;xlo,40,xup,40;xlo,1,xup,3;xlo,3,xup,3;xlo,5,xup,5;xlo,7,xup,7;xlo,9,xup,9;xlo,11,xup,11;xlo,13,xup,13;xlo,15,xup,15;xlo,17,xup,17;xlo,19,xup,19;xlo,21,xup,21;xlo,23,xup,23;xlo,25,xup,25;xlo,27,xup,27;xlo,29,xup,29;xlo,31,xup,31;xlo,33,xup,33;xlo,35,xup,35;xlo,37,xup,37;xlo,39,xup,39];

VV=vertcat(VV,borders);
%VV=vertcat(VV,borders,HORIZONTALS);
beginnings = horzcat(beginnings,bordersbeginnings');
endings=horzcat(endings,bordersendings');

%VV=JOINTS.X1Y1X2Y2;
thresold=input('Maximum length of rock bridges to fracture [m]? ');


[SurfArea,cycleList] = bridges_area(VV,thresold);
figure
%histpl=round(Nlines/2);
histogram(SurfArea)
mean(SurfArea)

%plot areas NOTE: THIS STEP IS TIME AND RAM CONSUMING
lengthplot=length(outx);

plotareas=figure(1);
%set(plotareas,'visible','off');
   plot(beginnings,endings,'color','black','linewidth',1);
   axis square
   axis([0,30,0,30])
   set(gca,'FontSize',14);
   hold on
for i=1:lengthplot
colfill=rand(1,3);
patch(outx{i,:}',outy{i,:}',colfill,'edgecolor','red','linewidth',1)

end
plot(beginnings,endings,'color','black','linewidth',1);
print(plotareas,'-djpeg','AREAPLOT')
toc
save ('SurfArea','SurfArea')
median(SurfArea)