% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% newlist=remove(list,entry);
%
function newlist = remove(list, entry);
newlist = list(find(list-entry*ones(size(list))));
