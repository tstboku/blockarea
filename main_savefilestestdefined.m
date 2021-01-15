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
JOINTS.SetA.amean=0;
JOINTS.SetA.astandard=5;
%random=0.5+13*rand(1,1);
%X=['Tracelength=',num2str(random)];
%disp(X)
JOINTS.SetA.lmean=2*xup;
random=3.3;
X=['Frequency=',num2str(random)];
disp(X)
JOINTS.SetA.onedimfreq=random; 
JOINTS.SetA.njoint=ceil(JOINTS.SetA.onedimfreq*xup);
JOINTS.SetA.angle=normrnd(JOINTS.SetA.amean,JOINTS.SetA.astandard,JOINTS.SetA.njoint,1);
% 2D joint point generation (Poisson process)
JOINTS.SetA.x1=rand(JOINTS.SetA.njoint,1);
JOINTS.SetA.length= (xup/2)*ones(JOINTS.SetA.njoint,1);


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

JOINTS.SetA.length= JOINTS.SetA.lmean*ones(JOINTS.SetA.njoint,1);
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

JOINTS.SetA.j1=find(JOINTS.SetA.x2<xlo);
JOINTS.SetA.dxnj=xlo-JOINTS.SetA.x1(JOINTS.SetA.j1,1);
JOINTS.SetA.dynj=JOINTS.SetA.dxnj.*tan(JOINTS.SetA.angle(JOINTS.SetA.j1,1)*pi/180);
JOINTS.SetA.dx(JOINTS.SetA.j1,1)=JOINTS.SetA.dxnj;
JOINTS.SetA.dy(JOINTS.SetA.j1,1)=JOINTS.SetA.dynj;
JOINTS.SetA.x2=JOINTS.SetA.x1+JOINTS.SetA.dx;
JOINTS.SetA.y2=JOINTS.SetA.y1+JOINTS.SetA.dy;

JOINTS.SetA.k1=find(JOINTS.SetA.y2>yup);
JOINTS.SetA.dynk=yup-JOINTS.SetA.y1(JOINTS.SetA.k1,1);
JOINTS.SetA.dxnk=JOINTS.SetA.dynk./tan(JOINTS.SetA.angle(JOINTS.SetA.k1,1)*pi/180);
JOINTS.SetA.dx(JOINTS.SetA.k1,1)=JOINTS.SetA.dxnk;
JOINTS.SetA.dy(JOINTS.SetA.k1,1)=JOINTS.SetA.dynk;
JOINTS.SetA.y2=JOINTS.SetA.y1+JOINTS.SetA.dy;
JOINTS.SetA.x2=JOINTS.SetA.x1+JOINTS.SetA.dx;

JOINTS.SetA.l1=find(JOINTS.SetA.y2<ylo);
JOINTS.SetA.dynl=ylo-JOINTS.SetA.y1(JOINTS.SetA.l1,1);
JOINTS.SetA.dxnl=JOINTS.SetA.dynl./tan(JOINTS.SetA.angle(JOINTS.SetA.l1,1)*pi/180);
JOINTS.SetA.dx(JOINTS.SetA.l1,1)=JOINTS.SetA.dxnl;
JOINTS.SetA.dy(JOINTS.SetA.l1,1)=JOINTS.SetA.dynl;
JOINTS.SetA.y2=JOINTS.SetA.y1+JOINTS.SetA.dy;
JOINTS.SetA.x2=JOINTS.SetA.x1+JOINTS.SetA.dx;

% cutting x3,y3
JOINTS.SetA.m1=find(JOINTS.SetA.x3>xup);
JOINTS.SetA.dxnm=JOINTS.SetA.x1(JOINTS.SetA.m1,1)-xup;
JOINTS.SetA.dynm=JOINTS.SetA.dxnm.*tan(JOINTS.SetA.angle(JOINTS.SetA.m1,1)*pi/180);
JOINTS.SetA.dx(JOINTS.SetA.m1,1)=JOINTS.SetA.dxnm;
JOINTS.SetA.dy(JOINTS.SetA.m1,1)=JOINTS.SetA.dynm;
JOINTS.SetA.x3=JOINTS.SetA.x1-JOINTS.SetA.dx;
JOINTS.SetA.y3=JOINTS.SetA.y1-JOINTS.SetA.dy;

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

JOINTS.SetA.p1=find(JOINTS.SetA.y3>yup);
JOINTS.SetA.dynp=JOINTS.SetA.y1(JOINTS.SetA.p1,1)-yup;
JOINTS.SetA.dxnp=JOINTS.SetA.dynp./tan(JOINTS.SetA.angle(JOINTS.SetA.p1,1)*pi/180);
JOINTS.SetA.dx(JOINTS.SetA.p1,1)=JOINTS.SetA.dxnp;
JOINTS.SetA.dy(JOINTS.SetA.p1,1)=JOINTS.SetA.dynp;
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



%startorient=90;
 %   random=-10+20*rand(1,1);  
 


JOINTS.SetB.astandard=5;
random=15;
X=['Tracelength=',num2str(random)];
disp(X)
JOINTS.SetB.lmean=random;
random=1.316761198;
X=['Frequency=',num2str(random)];
disp(X)
JOINTS.SetB.onedimfreq=random; 
JOINTS.SetB.njoint=ceil(JOINTS.SetB.onedimfreq*xup*(yup/JOINTS.SetB.lmean));
JOINTS.SetB.amean=90.99071697;
%JOINTS.SetB.angle=normrnd(JOINTS.SetB.amean,JOINTS.SetB.astandard,JOINTS.SetB.njoint,1);
% 2D joint point generation (Poisson process)
JOINTS.SetB.x1=xlo+(xup-xlo)*rand(JOINTS.SetB.njoint,1);
JOINTS.SetB.y1=ylo+(yup-ylo)*rand(JOINTS.SetB.njoint,1);
%JOINTS.SetB.angle=normrnd(JOINTS.SetB.amean,JOINTS.SetB.astandard,JOINTS.SetB.njoint,1);
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

JOINTS.SetB.j1=find(JOINTS.SetB.x2<xlo);
JOINTS.SetB.dxnj=xlo-JOINTS.SetB.x1(JOINTS.SetB.j1,1);
JOINTS.SetB.dynj=JOINTS.SetB.dxnj.*tan(JOINTS.SetB.angle(JOINTS.SetB.j1,1)*pi/180);
JOINTS.SetB.dx(JOINTS.SetB.j1,1)=JOINTS.SetB.dxnj;
JOINTS.SetB.dy(JOINTS.SetB.j1,1)=JOINTS.SetB.dynj;
JOINTS.SetB.x2=JOINTS.SetB.x1+JOINTS.SetB.dx;
JOINTS.SetB.y2=JOINTS.SetB.y1+JOINTS.SetB.dy;

JOINTS.SetB.k1=find(JOINTS.SetB.y2>yup);
JOINTS.SetB.dynk=yup-JOINTS.SetB.y1(JOINTS.SetB.k1,1);
JOINTS.SetB.dxnk=JOINTS.SetB.dynk./tan(JOINTS.SetB.angle(JOINTS.SetB.k1,1)*pi/180);
JOINTS.SetB.dx(JOINTS.SetB.k1,1)=JOINTS.SetB.dxnk;
JOINTS.SetB.dy(JOINTS.SetB.k1,1)=JOINTS.SetB.dynk;
JOINTS.SetB.y2=JOINTS.SetB.y1+JOINTS.SetB.dy;
JOINTS.SetB.x2=JOINTS.SetB.x1+JOINTS.SetB.dx;

JOINTS.SetB.l1=find(JOINTS.SetB.y2<ylo);
JOINTS.SetB.dynl=ylo-JOINTS.SetB.y1(JOINTS.SetB.l1,1);
JOINTS.SetB.dxnl=JOINTS.SetB.dynl./tan(JOINTS.SetB.angle(JOINTS.SetB.l1,1)*pi/180);
JOINTS.SetB.dx(JOINTS.SetB.l1,1)=JOINTS.SetB.dxnl;
JOINTS.SetB.dy(JOINTS.SetB.l1,1)=JOINTS.SetB.dynl;
JOINTS.SetB.y2=JOINTS.SetB.y1+JOINTS.SetB.dy;
JOINTS.SetB.x2=JOINTS.SetB.x1+JOINTS.SetB.dx;

% cutting x3,y3
JOINTS.SetB.m1=find(JOINTS.SetB.x3>xup);
JOINTS.SetB.dxnm=JOINTS.SetB.x1(JOINTS.SetB.m1,1)-xup;
JOINTS.SetB.dynm=JOINTS.SetB.dxnm.*tan(JOINTS.SetB.angle(JOINTS.SetB.m1,1)*pi/180);
JOINTS.SetB.dx(JOINTS.SetB.m1,1)=JOINTS.SetB.dxnm;
JOINTS.SetB.dy(JOINTS.SetB.m1,1)=JOINTS.SetB.dynm;
JOINTS.SetB.x3=JOINTS.SetB.x1-JOINTS.SetB.dx;
JOINTS.SetB.y3=JOINTS.SetB.y1-JOINTS.SetB.dy;

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

JOINTS.SetB.p1=find(JOINTS.SetB.y3>yup);
JOINTS.SetB.dynp=JOINTS.SetB.y1(JOINTS.SetB.p1,1)-yup;
JOINTS.SetB.dxnp=JOINTS.SetB.dynp./tan(JOINTS.SetB.angle(JOINTS.SetB.p1,1)*pi/180);
JOINTS.SetB.dx(JOINTS.SetB.p1,1)=JOINTS.SetB.dxnp;
JOINTS.SetB.dy(JOINTS.SetB.p1,1)=JOINTS.SetB.dynp;
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
%HORIZONTALS=[xlo,2,xup,2;xlo,4,xup,4;xlo,6,xup,6;xlo,8,xup,8;xlo,10,xup,10;xlo,12,xup,12;xlo,14,xup,14;xlo,16,xup,16;xlo,18,xup,18;xlo,20,xup,20;xlo,22,xup,22;xlo,24,xup,24;xlo,26,xup,26;xlo,28,xup,28;xlo,30,xup,30;xlo,32,xup,32;xlo,34,xup,34;xlo,36,xup,36;xlo,38,xup,38;xlo,40,xup,40;xlo,1,xup,3;xlo,3,xup,3;xlo,5,xup,5;xlo,7,xup,7;xlo,9,xup,9;xlo,11,xup,11;xlo,13,xup,13;xlo,15,xup,15;xlo,17,xup,17;xlo,19,xup,19;xlo,21,xup,21;xlo,23,xup,23;xlo,25,xup,25;xlo,27,xup,27;xlo,29,xup,29;xlo,31,xup,31;xlo,33,xup,33;xlo,35,xup,35;xlo,37,xup,37;xlo,39,xup,39];

VV=vertcat(VV,borders);
%VV=vertcat(VV,borders,HORIZONTALS);
beginnings = horzcat(beginnings,bordersbeginnings');
endings=horzcat(endings,bordersendings');

%VV=JOINTS.X1Y1X2Y2;
thresold=0.0;


%[SurfArea,cycleList] = bridges_area(VV,thresold);
[SurfArea,cycleList] = bridges_area(VV,thresold);


spacingA=1/JOINTS.SetA.onedimfreq;
spacingB=1/JOINTS.SetB.onedimfreq;
isbdwang10=spacingA*spacingB*0.479467001;
isbdwang20=spacingA*spacingB*0.795864071;
isbdwang30=spacingA*spacingB*1.133630478;
isbdwang40=spacingA*spacingB*1.508090213;
isbdwang50=spacingA*spacingB*1.942819438;
isbdwang60=spacingA*spacingB*2.511435611;
isbdwang70=spacingA*spacingB*3.252950031;
isbdwang80=spacingA*spacingB*4.310066626;
isbdwang90=spacingA*spacingB*6.171620416;
isbdwang100=spacingA*spacingB*11.48497618;
ISBD_Wang=[isbdwang10; isbdwang20; isbdwang30; isbdwang40; isbdwang50; isbdwang60; isbdwang70; isbdwang80; isbdwang90; isbdwang100];

sort_isbd_matlab=sort(datasample(SurfArea,5000));
cumsumisbd=cumsum(sort_isbd_matlab);
maxisbd=max(cumsumisbd);
n=numel(cumsumisbd);
for i=1:n
    yaxis_isbdmatlab(1,i)=(cumsumisbd(1,i)/maxisbd)*100;
end

val=10;
temp=abs(yaxis_isbdmatlab-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab(idx);
isbdmatlab10=sort_isbd_matlab(1,idx);

val=20;
temp=abs(yaxis_isbdmatlab-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab(idx);
isbdmatlab20=sort_isbd_matlab(1,idx);

val=30;
temp=abs(yaxis_isbdmatlab-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab(idx);
isbdmatlab30=sort_isbd_matlab(1,idx);

val=40;
temp=abs(yaxis_isbdmatlab-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab(idx);
isbdmatlab40=sort_isbd_matlab(1,idx);

val=50;
temp=abs(yaxis_isbdmatlab-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab(idx);
isbdmatlab50=sort_isbd_matlab(1,idx);

val=60;
temp=abs(yaxis_isbdmatlab-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab(idx);
isbdmatlab60=sort_isbd_matlab(1,idx);

val=70;
temp=abs(yaxis_isbdmatlab-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab(idx);
isbdmatlab70=sort_isbd_matlab(1,idx);

val=80;
temp=abs(yaxis_isbdmatlab-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab(idx);
isbdmatlab80=sort_isbd_matlab(1,idx);

val=90;
temp=abs(yaxis_isbdmatlab-val);
[idx idx]=min(temp);
closest=yaxis_isbdmatlab(idx);
isbdmatlab90=sort_isbd_matlab(1,idx);

isbdmatlab100=max(sort_isbd_matlab);
ISBD_NONPERSISTENT=[isbdmatlab10; isbdmatlab20; isbdmatlab30; isbdmatlab40; isbdmatlab50; isbdmatlab60; isbdmatlab70; isbdmatlab80; isbdmatlab90; isbdmatlab100];


isbdmult10=isbdmatlab10/isbdwang10;
isbdmult20=isbdmatlab20/isbdwang20;
isbdmult30=isbdmatlab30/isbdwang30;
isbdmult40=isbdmatlab40/isbdwang40;
isbdmult50=isbdmatlab50/isbdwang50;
isbdmult60=isbdmatlab60/isbdwang60;
isbdmult70=isbdmatlab70/isbdwang70;
isbdmult80=isbdmatlab80/isbdwang80;
isbdmult90=isbdmatlab90/isbdwang90;
isbdmult100=isbdmatlab100/isbdwang100;
ISBD_MULTIPLICATOR=[isbdmult10; isbdmult20; isbdmult30; isbdmult40; isbdmult50; isbdmult60; isbdmult70; isbdmult80; isbdmult90; isbdmult100];
    
SIZE_of_ANALYSIS=xup;
DIP_SetA=JOINTS.SetA.amean;
STD_of_DIP_SetA=JOINTS.SetA.astandard;
TRACELENGTH_SetA=JOINTS.SetA.lmean;
FREQUENCY_SetA=JOINTS.SetA.onedimfreq;
%DIP_SetB=mean(JOINTS.SetB.angle);
DIP_SetB=JOINTS.SetB.amean;
STD_of_DIP_SetB=JOINTS.SetB.astandard;
TRACELENGTH_SetB=JOINTS.SetB.lmean;
FREQUENCY_SetB=JOINTS.SetB.onedimfreq;

Z=isbdmult60;
ANGLE=abs(DIP_SetA-DIP_SetB);
X=sind(ANGLE)*JOINTS.SetA.lmean/JOINTS.SetB.onedimfreq;
Y=sind(ANGLE)*JOINTS.SetB.lmean/JOINTS.SetA.onedimfreq;

surfareatable=table(sort(SurfArea)');
isbdtable=table(ISBD_Wang,ISBD_NONPERSISTENT,ISBD_MULTIPLICATOR);
infotable=table(SIZE_of_ANALYSIS,DIP_SetA,STD_of_DIP_SetA,TRACELENGTH_SetA,FREQUENCY_SetA,DIP_SetB,STD_of_DIP_SetB,TRACELENGTH_SetB,FREQUENCY_SetB);
fittingtable=table(X,Y,Z);
filepath='D:/ISBDMATLAB/isbdinfo%s';
filexls='.xls';
%filename=sprintf(filepath,index,filexls);

filename2=[genvarname(['onepersonenon', num2str(datestr8601(now,'ymdHMS'))])];
filename=[filename2 '.xlsx'];
warning('off','MATLAB:xlswrite:AddSheet');
writetable(infotable,filename,'Sheet',1,'Range','A1');
writetable(isbdtable,filename,'Sheet',2,'Range','A1');
writetable(surfareatable,filename,'Sheet',3,'Range','A1');
writetable(fittingtable,filename,'Sheet',4,'Range','A1');
index=index+1;
clear beginnings borders bordersbeginnings berdersendings colfill cycleList endings fillings HORIZONTALS i JOINTS lengthplot  outx outy plotareas random setdescr start1  SurfArea thresold VV VV1 VV2 X
JOINTS.beginnings=[];
JOINTS.endings=[];
fillings=[];
end




















