function q=potential(L,b,alpha,x)
m=length(x);
q=zeros(m);

for k=1:m
	if abs(x)<L
		q(k)=1;
	else
		q(k)=b*exp(i*alpha);
	end
end

