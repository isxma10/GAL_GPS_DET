function [code_i,code_q]=L5codegen(svnum,nn,chip_rate,fs,code_delay)
%the code is represented at levels:-1 for bit=1;
%                                   1 for bit=0;
load L5I_codes.mat L5I_codes;
load L5Q_codes.mat L5Q_codes;
code1=L5I_codes(svnum,:);
code2=L5Q_codes(svnum,:);
ts=1/fs;           
% digitalize
code1_extend=[code1 code1 code1 code1 code1 code1 code1 code1 code1 code1];
code1_extend=[code1_extend code1_extend code1_extend code1_extend code1_extend];
code1_extend=[code1_extend code1_extend code1_extend code1_extend code1_extend];
code1_extend=[code1_extend code1_extend code1_extend code1_extend code1_extend];
code1_extend=[code1_extend code1_extend];
code2_extend=[code2 code2 code2 code2 code2 code2 code2 code2 code2 code2];
code2_extend=[code2_extend code2_extend code2_extend code2_extend code2_extend];
code2_extend=[code2_extend code2_extend code2_extend code2_extend code2_extend];
code2_extend=[code2_extend code2_extend code2_extend code2_extend code2_extend];
code2_extend=[code2_extend code2_extend];
b=[1:nn];
c=ceil(chip_rate*b*ts+10230-code_delay); 

code_i=code1_extend(c);
code_q=code2_extend(c);


