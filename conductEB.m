% Conduct permutation exchangable block.
% Code revised from hcp2blocks.m by Anderson M. Winkler
% http://brainder.org
% Reference:
% Winkler AM, Webster MA, Vidaurre D, Nichols TE, Smith SM.
% Multi-level block permutation. Neuroimage. 2015;123:253-68.

%load data
load('familyinfo.mat')

N = size(ID,1);

%sibtype
sibtype = zeros(N,1);
for n = 1:N,
    if famRela(n)<2,%single,sibling(not twin)=10
        sibtype(n) = 10;
    elseif famRela(n)==2,
        sibtype(n) = 100;%twin
    elseif famRela(n)==3,
        sibtype(n) = 1000;%triplet
    end
end

%famtype
Ufam = unique(famID);
famtype = zeros(N,1);
for f = 1:numel(Ufam),
    fidx = Ufam(f) == famID;
    famtype(fidx) = sum(sibtype(fidx));
end

%special subs, two pair of twins in one family
%no. 1607 and 3017 are twins
sibtype(1607)=101;
sibtype(3017)=101;

[~,idx] = sortrows([famID sibtype]);% ascending sort by famid(1) sibtype(2) age(3)
[~,idxback] = sort(idx);
sibtype = sibtype(idx);
famID   = famID(idx);
famtype = famtype(idx);

% Now make the blocks for each family
B = cell(numel(Ufam),1);
for f = 1:numel(Ufam),
    fidx = Ufam(f) == famID;
    ft = famtype(find(fidx,1));
    if any(ft == [10,20,200,400,3000]),
        B{f} = horzcat(famID(fidx),sibtype(fidx));
    else
        B{f} = horzcat(-famID(fidx),sibtype(fidx));
    end
end

% Concatenate all. Prepending the famtype ensures that the
% families of the same type can be shuffled whole-block. Also,
% add column with -1, for within-block at the outermost level
B = horzcat(-ones(N,1),famtype,cell2mat(B));

B = B(idxback,:);
            
save EB_abcd B ID;
%%
%subject No. column
B(:,5) = 1:1:3515;
[Pset, VG] = palm_quickperms([], B,5000);
save permOrder Pset VG ID;
            

