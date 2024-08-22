function ruo_trans2bin(X)

if(~exist("./data_bin/",'dir'))
        mkdir(char("./data_bin/"));
end
Fid = fopen("./data_bin/transmit_test.bin",'w');
fwrite(Fid,X,'int16');
% fclose('all');
fclose(Fid);
stop_addr = length(X)/704*176;
Fid2 = fopen("./data_bin/ddr_stop_addr_test.bin",'w');
fwrite(Fid2,stop_addr,'int32');
fclose('all');