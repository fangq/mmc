function hp=besselhprime(n,k,z)
hp=besselh(n-1,k,z)-n/z.*besselh(n,k,z);
