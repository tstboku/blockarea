function VVunicum = unicum(VV,Tolerance)

% delete the same and loops(points)

minx = min(min(VV(:,[1 3])));
maxx = max(max(VV(:,[1 3])));
miny = min(min(VV(:,[2 4])));
maxy = max(max(VV(:,[2 4])));

%Tolerance = 20*eps*max([maxy-miny maxx-minx abs([minx maxx miny maxy])]);
LogTol = ceil(log10(Tolerance));
VV = round(VV*10^(-LogTol))*10^LogTol;
VV = unique(VV,'rows','legacy');
[fr, ot] = ismember(VV,[VV(:,3:4) VV(:,1:2)],'rows','legacy');
ifr = find(ot>0);
Nfr = size(ifr,1);
if Nfr>0
    VVdelete = [];
    for i = 1:Nfr
        idelete = ot(ifr(i));
        if idelete>0
            VVdelete = [VVdelete; VV(idelete,:)];
            ot([ifr(i) ; idelete]) = 0;
        end;
    end;
    VV = setdiff(VV,VVdelete,'rows','legacy');
end;
VVunicum = VV;