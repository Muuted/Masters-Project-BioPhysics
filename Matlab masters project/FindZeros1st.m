%This function detects zeros crossings in the dataset (x,y) by linear
%interpolation between the cross points. Returns interpolated x-coordinate
%of zero crossing.
%
% FindZeros1st does not filter for tolerance 
%
function ZC = FindZeros1st(x,y)
zxidx=ZeroIndices(y); % get indices to zero crossings
ZC=x(zxidx); % do not interpolate. Just take smallest absolute(psi angle) before crossing zero
    function zci=ZeroIndices(v)
        zci1=find(v(:).*circshift(v(:), [-1 0]) <= 0); % Returns all zero-crossing indices by detecting sign change
        zci1=zci1(zci1~=numel(v) & zci1~=1); % removing endpoints from the zero crossing
        if length(zci1)>2 % only check for monotonous crossings if more than 2 crossings
            idm=abs(v(zci1-1))>abs(v(zci1)) & abs(v(zci1+2))>abs(v(zci1+1)); % return indices in zci1 where crossing is monotonous
            zci=zci1(idm); % return indices in input vector v with monotonous zero crossings.
        else
            zci=zci1;
        end
    end
end