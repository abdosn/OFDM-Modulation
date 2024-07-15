function recvoicedata( RxFileName , data_bin , quan_Levels ,nLevels , fs)


rec_im1 = reshape(data_bin,nLevels,numel(data_bin)/nLevels)'; 
rec_im2 = bin2dec(rec_im1); 
rec_im2(find(~rec_im2)) = 1;

recoverd = quan_Levels(rec_im2);

audiowrite(RxFileName,recoverd,fs);


end