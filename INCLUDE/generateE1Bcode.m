function E1Bcode = generateE1Bcode(PRN,BOC)

codes = load('codes_E1B.mat');
E1B=codes.codes_E1B(:,PRN);
chips=4092;
if BOC==1
    sub=ones(1,2*chips);
    E1Bdouble=ones(1,2*chips);
    for i=1:(2*chips)
        if rem(i,2) == 0
            sub(i)= -1;
        end %if
            E1Bdouble(i) = E1B(ceil(i/2));
            E1Bcode(i) = E1Bdouble(i)*sub(i);
    end %for
else 
   E1Bcode = E1B;
end %if

    