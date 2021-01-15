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
xup=50;
ylo=0;
yup=50;
%cd('D:\ISBDMATLAB');

JOINTS.beginnings=[];
JOINTS.endings=[];
% Number of joint sets

njointsets=2;
index=1;

for index=1:1000
startorient=0;
    random=-10+20*rand(1,1);  
JOINTS.SetA.amean=startorient+random;
JOINTS.SetA.astandard=5;
random=1+10*rand(1,1);
X=['Tracelength=',num2str(random)];
disp(X)
JOINTS.SetA.lmean=random;
random=0.5+3*rand(1,1);
X=['Frequency=',num2str(random)];
disp(X)
JOINTS.SetA.onedimfreq=random; 
JOINTS.SetA.njoint=ceil(JOINTS.SetA.onedimfreq*xup*(yup/JOINTS.SetA.lmean));
JOINTS.SetA.angle=normrnd(JOINTS.SetA.amean,JOINTS.SetA.astandard,JOINTS.SetA.njoint,1);
% 2D joint point generation (Poisson process)
JOINTS.SetA.x1=xlo+(xup-xlo)*rand(JOINTS.SetA.njoint,1);
JOINTS.SetA.y1=ylo+(yup-ylo)*rand(JOINTS.SetA.njoint,1);
JOINTS.SetA.angle=normrnd(JOINTS.SetA.amean,JOINTS.SetA.astandard,JOINTS.SetA.njoint,1);




    


% dip angle of joint
% Normal distribution
%JOINTS.SetA.amean=input('Mean dip angle of joints? ');
%JOINTS.SetA.astandard=input('Standard deviation of the dip angle? ');


% Simulation of the trace length distribution

% uniform distribution
% lmean=input('Trace length mean')
% random1=rand(njoint,1)
% length=(lmean*2).*random1

% neg. exponential distribution

JOINTS.SetA.length=exprnd(JOINTS.SetA.lmean,JOINTS.SetA.njoint,1);

% log. normal distribution
% lmean=input('Trace length mean')
% lstandard=input('standard deviation')
% length=lognrnd(lmean,lstandard,njoint,1)

% calculation of x2,y2
JOINTS.SetA.lh=JOINTS.SetA.length./2; %replace with lmean with length
JOINTS.SetA.dx=JOINTS.SetA.lh.*cos(JOINTS.SetA.angle*pi/180);
JOINTS.SetA.dy=JOINTS.SetA.lh.*sin(JOINTS.SetA.angle*pi/180);

JOINTS.SetA.x2=JOINTS.SetA.x1+JOINTS.SetA.dx;
JOINTS.SetA.y2=JOINTS.SetA.y1+JOINTS.SetA.dy;

JOINTS.SetA.x3=JOINTS.SetA.x1-JOINTS.SetA.dx;
JOINTS.SetA.y3=JOINTS.SetA.y1-JOINTS.SetA.dy;


% cutting x2,y2
JOINTS.SetA.i1=find(JOINTS.SetA.x2>xup);
JOINTS.SetA.dxni=xup-JOINTS.SetA.x1(JOINTS.SetA.i1,1);
JOINTS.SetA.dyni=JOINTS.SetA.dxni.*tan(JOINTS.SetA.angle(JOINTS.SetA.i1,1)*pi/180);
JOINTS.SetA.dx(JOINTS.SetA.i1,1)=JOINTS.SetA.dxni;
JOINTS.SetA.dy(JOINTS.SetA.i1,1)=JOINTS.SetA.dyni;
JOINTS.SetA.x2=JOINTS.SetA.x1+JOINTS.SetA.dx;
JOINTS.SetA.y2=JOINTS.SetA.y1+JOINTS.SetA.dy;



JOINTS.SetA.k1=find(JOINTS.SetA.y2>yup);
JOINTS.SetA.dynk=yup-JOINTS.SetA.y1(JOINTS.SetA.k1,1);
JOINTS.SetA.dxnk=JOINTS.SetA.dynk./tan(JOINTS.SetA.angle(JOINTS.SetA.k1,1)*pi/180);
JOINTS.SetA.dx(JOINTS.SetA.k1,1)=JOINTS.SetA.dxnk;
JOINTS.SetA.dy(JOINTS.SetA.k1,1)=JOINTS.SetA.dynk;
JOINTS.SetA.y2=JOINTS.SetA.y1+JOINTS.SetA.dy;
JOINTS.SetA.x2=JOINTS.SetA.x1+JOINTS.SetA.dx;



% cutting x3,y3


JOINTS.SetA.n1=find(JOINTS.SetA.x3<xlo);
JOINTS.SetA.dxnn=JOINTS.SetA.x1(JOINTS.SetA.n1,1)-xlo;
JOINTS.SetA.dynn=JOINTS.SetA.dxnn.*tan(JOINTS.SetA.angle(JOINTS.SetA.n1,1)*pi/180);
JOINTS.SetA.dx(JOINTS.SetA.n1,1)=JOINTS.SetA.dxnn;
JOINTS.SetA.dy(JOINTS.SetA.n1,1)=JOINTS.SetA.dynn;
JOINTS.SetA.x3=JOINTS.SetA.x1-JOINTS.SetA.dx;
JOINTS.SetA.y3=JOINTS.SetA.y1-JOINTS.SetA.dy;

JOINTS.SetA.o1=find(JOINTS.SetA.y3<ylo);
JOINTS.SetA.dyno=JOINTS.SetA.y1(JOINTS.SetA.o1,1)-ylo;
JOINTS.SetA.dxno=JOINTS.SetA.dyno./tan(JOINTS.SetA.angle(JOINTS.SetA.o1,1)*pi/180);
JOINTS.SetA.dx(JOINTS.SetA.o1,1)=JOINTS.SetA.dxno;
JOINTS.SetA.dy(JOINTS.SetA.o1,1)=JOINTS.SetA.dyno;
JOINTS.SetA.y3=JOINTS.SetA.y1-JOINTS.SetA.dy;
JOINTS.SetA.x3=JOINTS.SetA.x1-JOINTS.SetA.dx;



% new length
JOINTS.SetA.lnew=sqrt(JOINTS.SetA.dx.*JOINTS.SetA.dx+JOINTS.SetA.dy.*JOINTS.SetA.dy);
% total joint length
JOINTS.SetA.lsum=sum(JOINTS.SetA.lnew.*2);

% plot of points
% plot (x1,y1,'b+')
% hold on
%plot (x2,y2,'r*')

JOINTS.SetA.crackstring=zeros(JOINTS.SetA.njoint,1);
JOINTS.SetA.crackstring(:,1)=999999;
JOINTS.SetA.cracks=[JOINTS.SetA.crackstring,JOINTS.SetA.x3,JOINTS.SetA.y3,JOINTS.SetA.x2,JOINTS.SetA.y2];

% plot of fracture map
JOINTS.SetA.dxneu=JOINTS.SetA.x2-JOINTS.SetA.x3;
JOINTS.SetA.dyneu=JOINTS.SetA.y2-JOINTS.SetA.y3;
%quiver(x3,y3,dxneu,dyneu,0)
%plot(x3,y3)
%axis([xlo,xup,ylo,yup]);
%axis square;

JOINTS.SetA.beginnings=[JOINTS.SetA.x2 JOINTS.SetA.x3]';
JOINTS.SetA.endings=[JOINTS.SetA.y2 JOINTS.SetA.y3]';

%plot(JOINTS.SetA.beginnings, JOINTS.SetA.endings,'color','k');


JOINTS.beginnings=horzcat(JOINTS.beginnings,JOINTS.SetA.beginnings);
JOINTS.endings=horzcat(JOINTS.endings,JOINTS.SetA.endings);



startorient=90;
    random=-10+20*rand(1,1);  
 
JOINTS.SetB.amean=startorient+random;

JOINTS.SetB.astandard=5;
random=1+10*rand(1,1);
X=['Tracelength=',num2str(random)];
disp(X)
JOINTS.SetB.lmean=random;
random=0.5+3*rand(1,1);
X=['Frequency=',num2str(random)];
disp(X)
JOINTS.SetB.onedimfreq=random; 
JOINTS.SetB.njoint=ceil(JOINTS.SetB.onedimfreq*xup*(yup/JOINTS.SetB.lmean));
JOINTS.SetB.angle=normrnd(JOINTS.SetB.amean,JOINTS.SetB.astandard,JOINTS.SetB.njoint,1);
% 2D joint point generation (Poisson process)
JOINTS.SetB.x1=xlo+(xup-xlo)*rand(JOINTS.SetB.njoint,1);
JOINTS.SetB.y1=ylo+(yup-ylo)*rand(JOINTS.SetB.njoint,1);
JOINTS.SetB.angle=normrnd(JOINTS.SetB.amean,JOINTS.SetB.astandard,JOINTS.SetB.njoint,1);




    


% dip angle of joint
% Normal distribution
%JOINTS.SetB.amean=input('Mean dip angle of joints? ');
%JOINTS.SetB.astandard=input('Standard deviation of the dip angle? ');


% Simulation of the trace length distribution

% uniform distribution
% lmean=input('Trace length mean')
% random1=rand(njoint,1)
% length=(lmean*2).*random1

% neg. exponential distribution

JOINTS.SetB.length=exprnd(JOINTS.SetB.lmean,JOINTS.SetB.njoint,1);

% log. normal distribution
% lmean=input('Trace length mean')
% lstandard=input('standard deviation')
% length=lognrnd(lmean,lstandard,njoint,1)

% calculation of x2,y2
JOINTS.SetB.lh=JOINTS.SetB.length./2; %replace with lmean with length
JOINTS.SetB.dx=JOINTS.SetB.lh.*cos(JOINTS.SetB.angle*pi/180);
JOINTS.SetB.dy=JOINTS.SetB.lh.*sin(JOINTS.SetB.angle*pi/180);

JOINTS.SetB.x2=JOINTS.SetB.x1+JOINTS.SetB.dx;
JOINTS.SetB.y2=JOINTS.SetB.y1+JOINTS.SetB.dy;

JOINTS.SetB.x3=JOINTS.SetB.x1-JOINTS.SetB.dx;
JOINTS.SetB.y3=JOINTS.SetB.y1-JOINTS.SetB.dy;


% cutting x2,y2
JOINTS.SetB.i1=find(JOINTS.SetB.x2>xup);
JOINTS.SetB.dxni=xup-JOINTS.SetB.x1(JOINTS.SetB.i1,1);
JOINTS.SetB.dyni=JOINTS.SetB.dxni.*tan(JOINTS.SetB.angle(JOINTS.SetB.i1,1)*pi/180);
JOINTS.SetB.dx(JOINTS.SetB.i1,1)=JOINTS.SetB.dxni;
JOINTS.SetB.dy(JOINTS.SetB.i1,1)=JOINTS.SetB.dyni;
JOINTS.SetB.x2=JOINTS.SetB.x1+JOINTS.SetB.dx;
JOINTS.SetB.y2=JOINTS.SetB.y1+JOINTS.SetB.dy;



JOINTS.SetB.k1=find(JOINTS.SetB.y2>yup);
JOINTS.SetB.dynk=yup-JOINTS.SetB.y1(JOINTS.SetB.k1,1);
JOINTS.SetB.dxnk=JOINTS.SetB.dynk./tan(JOINTS.SetB.angle(JOINTS.SetB.k1,1)*pi/180);
JOINTS.SetB.dx(JOINTS.SetB.k1,1)=JOINTS.SetB.dxnk;
JOINTS.SetB.dy(JOINTS.SetB.k1,1)=JOINTS.SetB.dynk;
JOINTS.SetB.y2=JOINTS.SetB.y1+JOINTS.SetB.dy;
JOINTS.SetB.x2=JOINTS.SetB.x1+JOINTS.SetB.dx;



% cutting x3,y3


JOINTS.SetB.n1=find(JOINTS.SetB.x3<xlo);
JOINTS.SetB.dxnn=JOINTS.SetB.x1(JOINTS.SetB.n1,1)-xlo;
JOINTS.SetB.dynn=JOINTS.SetB.dxnn.*tan(JOINTS.SetB.angle(JOINTS.SetB.n1,1)*pi/180);
JOINTS.SetB.dx(JOINTS.SetB.n1,1)=JOINTS.SetB.dxnn;
JOINTS.SetB.dy(JOINTS.SetB.n1,1)=JOINTS.SetB.dynn;
JOINTS.SetB.x3=JOINTS.SetB.x1-JOINTS.SetB.dx;
JOINTS.SetB.y3=JOINTS.SetB.y1-JOINTS.SetB.dy;

JOINTS.SetB.o1=find(JOINTS.SetB.y3<ylo);
JOINTS.SetB.dyno=JOINTS.SetB.y1(JOINTS.SetB.o1,1)-ylo;
JOINTS.SetB.dxno=JOINTS.SetB.dyno./tan(JOINTS.SetB.angle(JOINTS.SetB.o1,1)*pi/180);
JOINTS.SetB.dx(JOINTS.SetB.o1,1)=JOINTS.SetB.dxno;
JOINTS.SetB.dy(JOINTS.SetB.o1,1)=JOINTS.SetB.dyno;
JOINTS.SetB.y3=JOINTS.SetB.y1-JOINTS.SetB.dy;
JOINTS.SetB.x3=JOINTS.SetB.x1-JOINTS.SetB.dx;



% new length
JOINTS.SetB.lnew=sqrt(JOINTS.SetB.dx.*JOINTS.SetB.dx+JOINTS.SetB.dy.*JOINTS.SetB.dy);
% total joint length
JOINTS.SetB.lsum=sum(JOINTS.SetB.lnew.*2);

% plot of points
% plot (x1,y1,'b+')
% hold on
%plot (x2,y2,'r*')

JOINTS.SetB.crackstring=zeros(JOINTS.SetB.njoint,1);
JOINTS.SetB.crackstring(:,1)=999999;
JOINTS.SetB.cracks=[JOINTS.SetB.crackstring,JOINTS.SetB.x3,JOINTS.SetB.y3,JOINTS.SetB.x2,JOINTS.SetB.y2];

% plot of fracture map
JOINTS.SetB.dxneu=JOINTS.SetB.x2-JOINTS.SetB.x3;
JOINTS.SetB.dyneu=JOINTS.SetB.y2-JOINTS.SetB.y3;
%quiver(x3,y3,dxneu,dyneu,0)
%plot(x3,y3)
%axis([xlo,xup,ylo,yup]);
%axis square;

JOINTS.SetB.beginnings=[JOINTS.SetB.x2 JOINTS.SetB.x3]';
JOINTS.SetB.endings=[JOINTS.SetB.y2 JOINTS.SetB.y3]';

%plot(JOINTS.SetB.beginnings, JOINTS.SetB.endings,'color','k');


JOINTS.beginnings=horzcat(JOINTS.beginnings,JOINTS.SetB.beginnings);
JOINTS.endings=horzcat(JOINTS.endings,JOINTS.SetB.endings);



JOINTS.X1Y1X2Y2=(vertcat(JOINTS.beginnings,JOINTS.endings))';


    VV1(:,1)=JOINTS.SetA.x2;
    VV1(:,2)=JOINTS.SetA.y2;
    VV1(:,3)=JOINTS.SetA.x3;
    VV1(:,4)=JOINTS.SetA.y3;
    VV2(:,1)=JOINTS.SetB.x2;
    VV2(:,2)=JOINTS.SetB.y2;
    VV2(:,3)=JOINTS.SetB.x3;
    VV2(:,4)=JOINTS.SetB.y3;
	VV=vertcat(VV1,VV2);


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



%[SurfArea,cycleList] = bridges_area(VV,thresold);
thresold=0.0;
[SurfArea,cycleList] = bridges_area(VV,thresold);
SurfArea0=SurfArea;
clear SurfArea;
thresold=0.05;
[SurfArea,cycleList] = bridges_area(VV,thresold);
SurfArea005=SurfArea;
clear SurfArea;
thresold=0.1;
[SurfArea,cycleList] = bridges_area(VV,thresold);
SurfArea01=SurfArea;
clear SurfArea;
thresold=0.2;
[SurfArea,cycleList] = bridges_area(VV,thresold);
SurfArea02=SurfArea;
clear SurfArea;
thresold=0.3;
[SurfArea,cycleList] = bridges_area(VV,thresold);
SurfArea03=SurfArea;
clear SurfArea;
thresold=0.4;
[SurfArea,cycleList] = bridges_area(VV,thresold);
SurfArea04=SurfArea;

%Table for threshold=0:
sort_isbd_matlab0=sort(datasample(SurfArea0,5000));
cumsumisbd0=cumsum(sort_isbd_matlab0);
maxisbd0=max(cumsumisbd0);
n=numel(cumsumisbd0);
for i=1:n
    yaxis_isbdmatlab0(1,i)=(cumsumisbd0(1,i)/maxisbd0)*100;
end

val=10;
temp=abs(yaxis_isbdmatlab0-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab0(idx);
isbdmatlab10_0=sort_isbd_matlab0(1,idx);

val=20;
temp=abs(yaxis_isbdmatlab0-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab0(idx);
isbdmatlab20_0=sort_isbd_matlab0(1,idx);

val=30;
temp=abs(yaxis_isbdmatlab0-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab0(idx);
isbdmatlab30_0=sort_isbd_matlab0(1,idx);

val=40;
temp=abs(yaxis_isbdmatlab0-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab0(idx);
isbdmatlab40_0=sort_isbd_matlab0(1,idx);

val=50;
temp=abs(yaxis_isbdmatlab0-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab0(idx);
isbdmatlab50_0=sort_isbd_matlab0(1,idx);

val=60;
temp=abs(yaxis_isbdmatlab0-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab0(idx);
isbdmatlab60_0=sort_isbd_matlab0(1,idx);

val=70;
temp=abs(yaxis_isbdmatlab0-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab0(idx);
isbdmatlab70_0=sort_isbd_matlab0(1,idx);

val=80;
temp=abs(yaxis_isbdmatlab0-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab0(idx);
isbdmatlab80_0=sort_isbd_matlab0(1,idx);

val=90;
temp=abs(yaxis_isbdmatlab0-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab0(idx);
isbdmatlab90_0=sort_isbd_matlab0(1,idx);

isbdmatlab100_0=max(sort_isbd_matlab0);
ISBD_0=[isbdmatlab10_0; isbdmatlab20_0; isbdmatlab30_0; isbdmatlab40_0; isbdmatlab50_0; isbdmatlab60_0; isbdmatlab70_0; isbdmatlab80_0; isbdmatlab90_0; isbdmatlab100_0];

%Table for threshold=0.05:
sort_isbd_matlab005=sort(datasample(SurfArea005,5000));
cumsumisbd005=cumsum(sort_isbd_matlab005);
maxisbd005=max(cumsumisbd005);
n=numel(cumsumisbd005);
for i=1:n
    yaxis_isbdmatlab005(1,i)=(cumsumisbd005(1,i)/maxisbd005)*100;
end

val=10;
temp=abs(yaxis_isbdmatlab005-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab005(idx);
isbdmatlab10_005=sort_isbd_matlab005(1,idx);

val=20;
temp=abs(yaxis_isbdmatlab005-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab005(idx);
isbdmatlab20_005=sort_isbd_matlab005(1,idx);

val=30;
temp=abs(yaxis_isbdmatlab005-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab005(idx);
isbdmatlab30_005=sort_isbd_matlab005(1,idx);

val=40;
temp=abs(yaxis_isbdmatlab005-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab005(idx);
isbdmatlab40_005=sort_isbd_matlab005(1,idx);

val=50;
temp=abs(yaxis_isbdmatlab005-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab005(idx);
isbdmatlab50_005=sort_isbd_matlab005(1,idx);

val=60;
temp=abs(yaxis_isbdmatlab005-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab005(idx);
isbdmatlab60_005=sort_isbd_matlab005(1,idx);

val=70;
temp=abs(yaxis_isbdmatlab005-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab005(idx);
isbdmatlab70_005=sort_isbd_matlab005(1,idx);

val=80;
temp=abs(yaxis_isbdmatlab005-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab005(idx);
isbdmatlab80_005=sort_isbd_matlab005(1,idx);

val=90;
temp=abs(yaxis_isbdmatlab005-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab005(idx);
isbdmatlab90_005=sort_isbd_matlab005(1,idx);

isbdmatlab100_005=max(SurfArea005);
ISBD_005=[isbdmatlab10_005; isbdmatlab20_005; isbdmatlab30_005; isbdmatlab40_005; isbdmatlab50_005; isbdmatlab60_005; isbdmatlab70_005; isbdmatlab80_005; isbdmatlab90_005; isbdmatlab100_005];

%Table for threshold=0.1:
sort_isbd_matlab01=sort(datasample(SurfArea01,5000));
cumsumisbd01=cumsum(sort_isbd_matlab01);
maxisbd01=max(cumsumisbd01);
n=numel(cumsumisbd01);
for i=1:n
    yaxis_isbdmatlab01(1,i)=(cumsumisbd01(1,i)/maxisbd01)*100;
end

val=10;
temp=abs(yaxis_isbdmatlab01-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab01(idx);
isbdmatlab10_01=sort_isbd_matlab01(1,idx);

val=20;
temp=abs(yaxis_isbdmatlab01-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab01(idx);
isbdmatlab20_01=sort_isbd_matlab01(1,idx);

val=30;
temp=abs(yaxis_isbdmatlab01-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab01(idx);
isbdmatlab30_01=sort_isbd_matlab01(1,idx);

val=40;
temp=abs(yaxis_isbdmatlab01-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab01(idx);
isbdmatlab40_01=sort_isbd_matlab01(1,idx);

val=50;
temp=abs(yaxis_isbdmatlab01-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab01(idx);
isbdmatlab50_01=sort_isbd_matlab01(1,idx);

val=60;
temp=abs(yaxis_isbdmatlab01-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab01(idx);
isbdmatlab60_01=sort_isbd_matlab01(1,idx);

val=70;
temp=abs(yaxis_isbdmatlab01-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab01(idx);
isbdmatlab70_01=sort_isbd_matlab01(1,idx);

val=80;
temp=abs(yaxis_isbdmatlab01-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab01(idx);
isbdmatlab80_01=sort_isbd_matlab01(1,idx);

val=90;
temp=abs(yaxis_isbdmatlab01-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab01(idx);
isbdmatlab90_01=sort_isbd_matlab01(1,idx);

isbdmatlab100_01=max(SurfArea01);
ISBD_01=[isbdmatlab10_01; isbdmatlab20_01; isbdmatlab30_01; isbdmatlab40_01; isbdmatlab50_01; isbdmatlab60_01; isbdmatlab70_01; isbdmatlab80_01; isbdmatlab90_01; isbdmatlab100_01];

%Table for threshold=0.2:
sort_isbd_matlab02=sort(datasample(SurfArea02,5000));
cumsumisbd02=cumsum(sort_isbd_matlab02);
maxisbd02=max(cumsumisbd02);
n=numel(cumsumisbd02);
for i=1:n
    yaxis_isbdmatlab02(1,i)=(cumsumisbd02(1,i)/maxisbd02)*100;
end

val=10;
temp=abs(yaxis_isbdmatlab02-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab02(idx);
isbdmatlab10_02=sort_isbd_matlab02(1,idx);

val=20;
temp=abs(yaxis_isbdmatlab02-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab02(idx);
isbdmatlab20_02=sort_isbd_matlab02(1,idx);

val=30;
temp=abs(yaxis_isbdmatlab02-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab02(idx);
isbdmatlab30_02=sort_isbd_matlab02(1,idx);

val=40;
temp=abs(yaxis_isbdmatlab02-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab02(idx);
isbdmatlab40_02=sort_isbd_matlab02(1,idx);

val=50;
temp=abs(yaxis_isbdmatlab02-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab02(idx);
isbdmatlab50_02=sort_isbd_matlab02(1,idx);

val=60;
temp=abs(yaxis_isbdmatlab02-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab02(idx);
isbdmatlab60_02=sort_isbd_matlab02(1,idx);

val=70;
temp=abs(yaxis_isbdmatlab02-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab02(idx);
isbdmatlab70_02=sort_isbd_matlab02(1,idx);

val=80;
temp=abs(yaxis_isbdmatlab02-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab02(idx);
isbdmatlab80_02=sort_isbd_matlab02(1,idx);

val=90;
temp=abs(yaxis_isbdmatlab02-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab02(idx);
isbdmatlab90_02=sort_isbd_matlab02(1,idx);

isbdmatlab100_02=max(SurfArea02);
ISBD_02=[isbdmatlab10_02; isbdmatlab20_02; isbdmatlab30_02; isbdmatlab40_02; isbdmatlab50_02; isbdmatlab60_02; isbdmatlab70_02; isbdmatlab80_02; isbdmatlab90_02; isbdmatlab100_02];

%Table for threshold=0.3:
sort_isbd_matlab03=sort(datasample(SurfArea03,5000));
cumsumisbd03=cumsum(sort_isbd_matlab03);
maxisbd03=max(cumsumisbd03);
n=numel(cumsumisbd03);
for i=1:n
    yaxis_isbdmatlab03(1,i)=(cumsumisbd03(1,i)/maxisbd03)*100;
end

val=10;
temp=abs(yaxis_isbdmatlab03-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab03(idx);
isbdmatlab10_03=sort_isbd_matlab03(1,idx);

val=20;
temp=abs(yaxis_isbdmatlab03-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab03(idx);
isbdmatlab20_03=sort_isbd_matlab03(1,idx);

val=30;
temp=abs(yaxis_isbdmatlab03-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab03(idx);
isbdmatlab30_03=sort_isbd_matlab03(1,idx);

val=40;
temp=abs(yaxis_isbdmatlab03-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab03(idx);
isbdmatlab40_03=sort_isbd_matlab03(1,idx);

val=50;
temp=abs(yaxis_isbdmatlab03-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab03(idx);
isbdmatlab50_03=sort_isbd_matlab03(1,idx);

val=60;
temp=abs(yaxis_isbdmatlab03-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab03(idx);
isbdmatlab60_03=sort_isbd_matlab03(1,idx);

val=70;
temp=abs(yaxis_isbdmatlab03-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab03(idx);
isbdmatlab70_03=sort_isbd_matlab03(1,idx);

val=80;
temp=abs(yaxis_isbdmatlab03-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab03(idx);
isbdmatlab80_03=sort_isbd_matlab03(1,idx);

val=90;
temp=abs(yaxis_isbdmatlab03-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab03(idx);
isbdmatlab90_03=sort_isbd_matlab03(1,idx);

isbdmatlab100_03=max(SurfArea03);
ISBD_03=[isbdmatlab10_03; isbdmatlab20_03; isbdmatlab30_03; isbdmatlab40_03; isbdmatlab50_03; isbdmatlab60_03; isbdmatlab70_03; isbdmatlab80_03; isbdmatlab90_03; isbdmatlab100_03];

%Table for threshold=0.4:
sort_isbd_matlab04=sort(datasample(SurfArea04,5000));
cumsumisbd04=cumsum(sort_isbd_matlab04);
maxisbd04=max(cumsumisbd04);
n=numel(cumsumisbd04);
for i=1:n
    yaxis_isbdmatlab04(1,i)=(cumsumisbd04(1,i)/maxisbd04)*100;
end

val=10;
temp=abs(yaxis_isbdmatlab04-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab04(idx);
isbdmatlab10_04=sort_isbd_matlab04(1,idx);

val=20;
temp=abs(yaxis_isbdmatlab04-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab04(idx);
isbdmatlab20_04=sort_isbd_matlab04(1,idx);

val=30;
temp=abs(yaxis_isbdmatlab04-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab04(idx);
isbdmatlab30_04=sort_isbd_matlab04(1,idx);

val=40;
temp=abs(yaxis_isbdmatlab04-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab04(idx);
isbdmatlab40_04=sort_isbd_matlab04(1,idx);

val=50;
temp=abs(yaxis_isbdmatlab04-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab04(idx);
isbdmatlab50_04=sort_isbd_matlab04(1,idx);

val=60;
temp=abs(yaxis_isbdmatlab04-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab04(idx);
isbdmatlab60_04=sort_isbd_matlab04(1,idx);

val=70;
temp=abs(yaxis_isbdmatlab04-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab04(idx);
isbdmatlab70_04=sort_isbd_matlab04(1,idx);

val=80;
temp=abs(yaxis_isbdmatlab04-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab04(idx);
isbdmatlab80_04=sort_isbd_matlab04(1,idx);

val=90;
temp=abs(yaxis_isbdmatlab04-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab04(idx);
isbdmatlab90_04=sort_isbd_matlab04(1,idx);

isbdmatlab100_04=max(SurfArea04);
ISBD_04=[isbdmatlab10_04; isbdmatlab20_04; isbdmatlab30_04; isbdmatlab40_04; isbdmatlab50_04; isbdmatlab60_04; isbdmatlab70_04; isbdmatlab80_04; isbdmatlab90_04; isbdmatlab100_04];

    
SIZE_of_ANALYSIS=xup;
DIP_SetA=JOINTS.SetA.amean;
STD_of_DIP_SetA=JOINTS.SetA.astandard;
TRACELENGTH_SetA=JOINTS.SetA.lmean;
FREQUENCY_SetA=JOINTS.SetA.onedimfreq;
DIP_SetB=JOINTS.SetB.amean;
STD_of_DIP_SetB=JOINTS.SetB.astandard;
TRACELENGTH_SetB=JOINTS.SetB.lmean;
FREQUENCY_SetB=JOINTS.SetB.onedimfreq;

diff005=isbdmatlab50_005/isbdmatlab50_0;
diff01=isbdmatlab50_01/isbdmatlab50_0;
diff02=isbdmatlab50_02/isbdmatlab50_0;
diff03=isbdmatlab50_03/isbdmatlab50_0;
diff04=isbdmatlab50_04/isbdmatlab50_0;

surfarea0=table(sort(SurfArea0)');
surfarea005=table(sort(SurfArea005)');
surfarea01=table(sort(SurfArea01)');
surfarea02=table(sort(SurfArea02)');
surfarea03=table(sort(SurfArea03)');
surfarea04=table(sort(SurfArea04)');
isbdtable=table(ISBD_0,ISBD_005,ISBD_01,ISBD_02,ISBD_03,ISBD_04);
differencetable=table(diff005,diff01,diff02,diff03,diff04);
infotable=table(SIZE_of_ANALYSIS,DIP_SetA,STD_of_DIP_SetA,TRACELENGTH_SetA,FREQUENCY_SetA,DIP_SetB,STD_of_DIP_SetB,TRACELENGTH_SetB,FREQUENCY_SetB);

filepath='D:/ISBDMATLAB/isbdinfo%s';
filexls='.xls';
%filename=sprintf(filepath,index,filexls);

filename2=[genvarname(['isbdinfo', num2str(datestr8601(now,'ymdHMS'))])];
filename=[filename2 '.xlsx'];
warning('off','MATLAB:xlswrite:AddSheet');
writetable(infotable,filename,'Sheet',1,'Range','A1');
writetable(isbdtable,filename,'Sheet',2,'Range','A1');
writetable(surfarea0,filename,'Sheet',3,'Range','A1');
writetable(surfarea005,filename,'Sheet',4,'Range','A1');
writetable(surfarea01,filename,'Sheet',5,'Range','A1');
writetable(surfarea02,filename,'Sheet',6,'Range','A1');
writetable(surfarea03,filename,'Sheet',7,'Range','A1');
writetable(surfarea04,filename,'Sheet',8,'Range','A1');
writetable(differencetable,filename,'Sheet',9,'Range','A1');
index=index+1;
clear beginnings borders bordersbeginnings bordersendings colfill cycleList endings fillings HORIZONTALS i JOINTS lengthplot  outx outy plotareas random setdescr start1  SurfArea thresold VV VV1 VV2 X
JOINTS.beginnings=[];
JOINTS.endings=[];
fillings=[];
end













