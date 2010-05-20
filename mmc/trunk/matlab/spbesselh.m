function hn=spbesselh(n,k,z)
if(k==1)
    hn=spbesselj(n,z)+i*spbessely(n,z);
elseif(k==2)
    hn=spbesselj(n,z)-i*spbessely(n,z);
else
    error('wrong value for the second parameter');
end
