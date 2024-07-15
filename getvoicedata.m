function [ser_data_bin , quan_Levels ,fs] = getvoicedata(FileName , nLevels)

[y, fs] = audioread(FileName);


[quan_Levels , x2 , ~ ] = cuantizar(y,1,(2^nLevels)-1);

bit1 = (dec2bin(x2));
bit2 = bit1.';
ser_data_bin = bit2(:).';
%{
rec_im1 = reshape(data_bin,nLevels,numel(data_bin)/nLevels)'; 
rec_im2 = bin2dec(rec_im1); 
rec_im3 = reshape(rec_im2,size(y));


recoverd = quan_Levels(rec_im3);

audiowrite('rec.wav',recoverd,fs);
%}
end