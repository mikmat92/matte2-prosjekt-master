%lagA Lager koeffisientmatrisen A for en Euler-Bernoullibjelke
%	Input:
%		n: antall likninger, dvs hvor mange segmenter vi deler bjelken vår
%			i, slik at h = L/n
function A = lagA( n )
	A = spdiags(repmat([1 -4 6 -4 1], n, 1),(-2:2),n,n);
	A(1,1:4)= [16 -9 8/3 -1/4];
	A(n-1,n-3:n)= [16 -60 72 -28]/17; 
	A(n ,n-3:n)= [-12 96 -156 72]/17;
end