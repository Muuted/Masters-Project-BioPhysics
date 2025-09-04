% 2nd zero detection function 
%
% FindZeros2nd filters for tolerance 
%
function [zc,idz] = FindZeros2nd(x,y)
idz=ZeroIndices(y); % get indices to zero crossings
%zc=x(idz); % do not interpolate. Just take smallest absolute(psi angle) before crossing zero
%
zc=zeros(1,numel(idz)); % initialize ZC
for k1 = 1:numel(idz) % procedure to interpolate
    idxrng = max([1 idz(k1)-1]):min([idz(k1)+1 numel(y)]);
    xrng = x(idxrng);
    yrng = y(idxrng);
    zc(k1) = interp1( yrng(:), xrng(:), 0, 'linear', 'extrap' );
end
%
    function zi=ZeroIndices(v)
        zi1=find(v(:).*circshift(v(:), [-1 0]) <= 0); % Find lower bound index to zero-crossing by detecting sign change
        zi1=zi1(zi1~=numel(v) & zi1~=1); % removing endpoints from the zero crossing
        if length(zi1)>2 % check for monotonous crossings if more than 2 crossings
            idm=abs(v(zi1-1))>abs(v(zi1)) & abs(v(zi1+2))>abs(v(zi1+1)); % return indices idm in zi1 where crossing is monotonous
            zi1=zi1(idm); % return indices in input vector v with monotonous zero crossings.
        end
        % filter zero crossing for tolerance
        tolerance=0.2; % tolerance in deviation of kG from target value
        idp=abs(v(zi1))<tolerance; % require better than tolerance
        zi1=zi1(idp);
        zi=zi1;
    end
end