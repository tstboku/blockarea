function VVunicum = unicum2(VV,Tolerance)

% delete the same and loops(points)

minx = min(min(VV(:,1)));
maxx = max(max(VV(:,1)));
miny = min(min(VV(:,2)));
maxy = max(max(VV(:,2)));

%Tolerance = 20*eps*max([maxy-miny maxx-minx abs([minx maxx miny maxy])]);
LogTol = ceil(log10(Tolerance));
VV = round(VV*10^(-LogTol))*10^LogTol;
VV = unique(VV,'rows','legacy');
VVunicum = VV;