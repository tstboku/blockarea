function [E,V,Nnew] = delete_trees(E,V)
% delete trees
iaa = [1];
while ~isempty(iaa)
    NV = size(V,1);
    [aa] = hist(E,1:NV);
    iaa = find(sum(aa,2)==1); 
    E = E(sum(ismember(E,iaa,'legacy'),2)==0,:);
    VE = unique([E(:,1); E(:,2)],'rows','legacy');
    [~,VeV] = ismember(1:NV,VE,'legacy');
    E = VeV(E);
    V = V(VE,:);
    Nnew = size(E,1);
end;
if Nnew == 1    
    E = [];V = []; Nnew = 0;
end;