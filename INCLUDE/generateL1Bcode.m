function L1Bcode = generateL1Bcode(chips)
% generateL1Bcode.m generates one of the PRN code transmitted by the GIOVE
% satellite on L1 B
%
fullCode= rot90(load('L1-B-1.dat'));
%--- Generate a defined number of chips of the L1B code -----
for i=1:chips    
    if fullCode(i)== 1
        L1Bcode(i) = -1;
    else
        L1Bcode(i) = 1;
    end %if
end %for

