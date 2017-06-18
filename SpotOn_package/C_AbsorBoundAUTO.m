function y=C_AbsorBoundAUTO(z, CurrTime, D, halfZ)
%This is a corrected version of equation 16 in Kues and Kubitscheck, Single
%Molecules, 2002 and a corrected version of equation Suppl 5.7 in Mazza et
%al, Nucleic Acids Research, 2012. Both equations are wrong, but they are
%wrong in different ways and the one below is correct. 
%Moreover, this implementation automatically stops the sum when the error
%is negligble. 

WhenToStop = 1e-10;
f = inf;
n=0; %iterator
h=1;
while abs(f) > WhenToStop
    f = (-1)^n * ( erfc( ((2*n+1)*halfZ-z)/sqrt(4*D*CurrTime)) +  erfc( ((2*n+1)*halfZ+z)/sqrt(4*D*CurrTime))  );
    h = h - f;
    n = n + 1;
end
y=h;    

end

