function [SurfArea,cycleList] = bridges_area(VV,thresold)
global fillings out2
polysx=[];
polysy=[];
sepa=NaN;
beginnings=[VV(:,1),VV(:,3)]';
endings=[VV(:,2),VV(:,4)]';





minx = min(min(VV(:,[1 3])));
maxx = max(max(VV(:,[1 3])));
miny = min(min(VV(:,[2 4])));
maxy = max(max(VV(:,[2 4])));

Tolerance = 200*eps*max([maxy-miny maxx-minx abs([minx maxx miny maxy])]);

VV = unicum(VV,Tolerance);
Nlines = size(VV,1);

%%

% intersection

VVbridge = [];
countBridge = 1;
VVnew = [];
for i = 1:Nlines  
    count_i = 1;
    VVi = [];
    p1 = VV(i,1:2); 
    p2 = VV(i,3:4);  
    out = lineSegmentIntersect(VV,VV(i,:));
    for j = 1:Nlines
        if i == j
            continue;
        end;
        p3 = VV(j,1:2); 
        p4 = VV(j,3:4);
        if out.intAdjacencyMatrix(j)
            % intersection
            Pint = [out.intMatrixX(j) out.intMatrixY(j)];
            VVi(count_i,:) = [Pint];
            count_i = count_i + 1;
        else
            % bridge
            [distance1,~,Vi1,Vj1,seg1] = DistBetween2Segment(p1,p2,p3,p4);
            [distance2,~,Vi2,Vj2,seg2] = DistBetween2Segment(p2,p1,p3,p4);
            if distance1<=thresold
                if seg1 == 2 || seg1 == 3
                    VVi(count_i,:) = [Vi1];
                    VVbridge(countBridge,:) = [Vi1,Vj1];
                    count_i = count_i + 1;
                    countBridge = countBridge + 1;
                end;
            end;
            if sum(Vi1-Vi2)~=0 & distance2<=thresold
                if seg2 == 2 || seg2 == 3
                    VVi(count_i,:) = [Vi2];
                    VVbridge(countBridge,:) = [Vi2,Vj2];
                    count_i = count_i + 1;
                    countBridge = countBridge + 1;
                end;
            end;
        end;        
    end;
    % adding all intersection and bridges
    if ~isempty(VVi)
        VVi = unique([p1; VVi; p2],'rows','legacy');
        iuniques = find(abs(sum(diff(VVi).^2,2))>Tolerance);
        VVi = [VVi(iuniques,:); VVi(end,:)];
        VVnew = [VVnew
            VVi(2:end,:) VVi(1:end-1,:)];
    else
        VVnew = [VVnew
            VV(i,:)];
    end;
end;

% adding bridges
VVV = [VVnew
        VVbridge];
VVV = unicum(VVV,Tolerance);

% intersection with bridges and between bridges
% intersection

VVnew = [];
for i = 1:size(VVV,1)  
    out = lineSegmentIntersect(VVV,VVV(i,:));
    J = find(out.intAdjacencyMatrix);
    if ~isempty(J)
        % intersections
        Pint = [out.intMatrixX(J) out.intMatrixY(J)];
        % adding all intersections
        p1 = VVV(i,1:2); 
        p2 = VVV(i,3:4); 
        VVi = unique([p1; Pint; p2],'rows','legacy');
        iuniques = find(abs(sum(diff(VVi).^2,2))>Tolerance);
        VVi = [VVi(iuniques,:); VVi(end,:)];
        VVnew = [VVnew
            VVi(2:end,:) VVi(1:end-1,:)];
    else
        VVnew = [VVnew
            VVV(i,:)];
    end;
end;

VVnew = unicum(VVnew,Tolerance);

% find all vertex
V = unique([VVnew(:,1:2); VVnew(:,3:4)],'rows','legacy');
V = unicum2(V,Tolerance);

% convert task to the graph task
Nnew = size(VVnew,1);
E = zeros(Nnew,2);

for i = 1:Nnew
    E(i,1) = find(V(:,1) == VVnew(i,1) & V(:,2) == VVnew(i,2));
    E(i,2) = find(V(:,1) == VVnew(i,3) & V(:,2) == VVnew(i,4));
end;

% delete trees
[E,V,Nnew] = delete_trees(E,V);

SurfArea = [];
cycleList = {};
icycle = 1;

%% loop over all boundaries
while ~isempty(E)
    % find boundary
    minX = min(V(:,1));
    iX0 = find(V(:,1) == minX,1);
    Pout = [iX0];
    iXprev = [-1];
    iXcurr = iX0;
    ibegin = find(Pout==iXcurr,1);
    angle_prev = atan2(1,0);     %pi/2
    % loop over
    while (size(Pout,2)==1) || (isempty(ibegin) && size(Pout,2)>1) 
        ineir = find(E==iXcurr);
        iXneir = E(mod(ineir+Nnew-1,2*Nnew)+1);
        iXneir = setdiff(iXneir,iXprev,'legacy');
        angle_next = atan2(V(iXneir,2)-V(iXcurr,2),V(iXneir,1)-V(iXcurr,1));
        diff_angle = mod(angle_next-angle_prev+pi,2*pi);
        diff_angle = round(diff_angle*1e13)*1e-13;
        [~,iiXneir] = max(diff_angle);
        if size(iiXneir,2)==1
            iXprev = iXcurr;
            iXcurr = iXneir(iiXneir);
            ibegin = find(Pout==iXcurr,1);
            Pout = [Pout iXcurr];
            angle_prev  = angle_next(iiXneir);
        else
            % should be nearst between with the same angles
            iiiXneir = min((V(iXneir(iiXneir),2)-V(iXcurr,2)).^2 + (V(iXneir(iiXneir),1)-V(iXcurr,1)).^2);
            iXprev = iXcurr;
            iXcurr = iXneir(iiXneir(iiiXneir));
            ibegin = find(Pout==iXcurr,1);
            Pout = [Pout iXcurr];
            angle_prev  = angle_next(iiXneir(iiiXneir));
        end;
    end;
    Pout = Pout(ibegin:end);
    
    Vnew = [];
    Enew = [];
    icycle = 1;    
    PPout{1} = Pout;
    Kout = 1;
    while ~isempty(PPout)
        Pout = PPout{1};
        Pin = Pout(1:2);
        iXprev = Pout(1);
        iXcurr = Pout(2);
        angle_prev = atan2(V(iXcurr,2)-V(iXprev,2),V(iXcurr,1)-V(iXprev,1));
        % loop over 
        while (iXcurr~=Pin(1) & size(Pin,2)>1) || size(Pin,2)==1
            ineir = find(E==iXcurr);
            iXneir = E(mod(ineir+Nnew-1,2*Nnew)+1);
            iXneir = setdiff(iXneir,iXprev,'legacy');
            angle_next = atan2(V(iXneir,2)-V(iXcurr,2),V(iXneir,1)-V(iXcurr,1));
            diff_angle = mod(angle_next-angle_prev+pi,2*pi);
            diff_angle = round(diff_angle*1e13)*1e-13;
            [~,iiXneir] = min(diff_angle);
            if size(iiXneir,2)==1
                iXprev = iXcurr;
                iXcurr = iXneir(iiXneir);
                Pin = [Pin iXcurr];
                angle_prev  = angle_next(iiXneir);
            else
                % should be nearst between with the same angles
                iiiXneir = min((V(iXneir(iiXneir),2)-V(iXcurr,2)).^2 + (V(iXneir(iiXneir),1)-V(iXcurr,1)).^2);
                iXprev = iXcurr;
                iXcurr = iXneir(iiXneir(iiiXneir));
                Pin = [Pin iXcurr];
                angle_prev  = angle_next(iiXneir(iiiXneir));
            end;
        end;
        cycleList{icycle} = Pin;
        SurfArea(icycle) = polyarea(V(Pin,1),V(Pin,2));   
       
           polysx=vertcat(polysx,V(Pin,1),sepa);
    polysy=vertcat(polysy,V(Pin,2),sepa);
        
        
        icycle = icycle + 1;
        Vnew = union(Vnew,Pin,'legacy');
        if isempty(Enew)
            Enew = [Pin(1:end-1)' Pin(2:end)'];
        else
            Enew = union(Enew,[Pin(1:end-1)' Pin(2:end)'],'rows','legacy');
        end;        
        
        
        [Pd,iPd] = ismember(Pout,Pin,'legacy');
        [Pd2,iPd2] = ismember(Pin,Pout,'legacy');
        iPin = find(Pd==0);
        if isempty(iPin)
            Pout = [];
            PPout = PPout(2:end);
            Kout = Kout - 1;
        else
            ibegin = find(diff(Pd) == -1);
            iend = find(diff(Pd) == 1);
            K = size(ibegin,2);
            k = 1;
            P = [ Pout(ibegin(k):iend(k)+1)   fliplr(Pin(iPd(ibegin(k))+1:iPd(iend(k)+1)-1))   Pout(ibegin(k))];
            PPout{1} = P;
            for k = 2:K
                P = [ Pout(ibegin(k):iend(k)+1)   fliplr(Pin(iPd(ibegin(k))+1:iPd(iend(k)+1)-1))   Pout(ibegin(k))];
                Kout = Kout + 1;
                PPout{Kout} = P;
            end;
        end;
    end;
    E = setdiff( setdiff(E,fliplr(Enew),'rows','legacy'),Enew,'rows','legacy');
    % delete trees
    [E,V,Nnew] = delete_trees(E,V);
    
    
    
    %PLOTTING RESULTS:
a=polysx';
if ~isnan(a(end))
  a(end+1)=nan;
end
idx2=find(isnan(a));
idx1=[1 idx2(1:end-1)+1];
n=numel(idx1);
outx=cell(n,1);
for k=1:n
  outx{k}=a(idx1(k):idx2(k)-1);
end
assignin('base', 'outx', outx);

b=polysy';
if ~isnan(b(end))
  b(end+1)=nan;
end
idx4=find(isnan(b));
idx3=[1 idx4(1:end-1)+1];
n=numel(idx3);
outy=cell(n,1);
for k=1:n
  outy{k}=b(idx3(k):idx4(k)-1);
end
assignin('base', 'outy', outy);
end;