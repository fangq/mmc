function hn=spbesselhprime(n,k,z)
if(k==1)
    hn=spbesseljprime(n,z)+i*spbesselyprime(n,z);
elseif(k==2)
    hn=spbesseljprime(n,z)-i*spbesselyprime(n,z);
else
    error('wrong value for the second parameter');
end