input_filename='PRN1_m130dBm_CH1.bin';
output_filename=[input_filename(1:end-4) '_int8.dat'];

offset = 2^10; 
fid_input = fopen(input_filename,'r');
fid_output = fopen(output_filename,'w');

fseek(fid_input,offset,'bof');

samplingFreq = 53e6;
blksize_sec = 0.1;
blksize = samplingFreq*blksize_sec;
samplesRead = blksize;
dataOut=zeros(1,blksize);
while (samplesRead == blksize)
    [data, samplesRead] = fread(fid_input,blksize,'ubit2');

    for i = 1:length(data)
        switch(data(i))
            case 2 
                dataOut(i) = 3;
            case 0 
                dataOut(i) = 1;
            case 1
                dataOut(i) = -1;
            case 3
                dataOut(i) = -3;        
        end 
    end
    % write the output data to a file
    fwrite(fid_output, dataOut, 'int8');
    
end
disp('End of file');
fclose(fid_input);
fclose(fid_output);

