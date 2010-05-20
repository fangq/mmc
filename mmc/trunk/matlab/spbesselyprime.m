function yp=spbesselyprime(n,z)
yp=besselyprime(n+1/2,z).*sqrt(pi/(2*z))-sqrt(pi/2)*bessely(n+1/2,z)./(2*z.*sqrt(z));
