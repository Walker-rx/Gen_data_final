% clear;
close all;

channel_choice = 1; 
channel_choice_inf = 3;
dir_up = "./data_set_final/";
test = [1 0 1 0 1 1 1];
% bias = input('bias(mA): ')

% dp821A = visa('ni','USB0::0x1AB1::0x0E11::DP8E163250125::INSTR');
% % dp821A = visa('ni','ASRL3:INSTR');
% fopen(dp821A);
% % fprintf(dp821A,'*IDN?'); 
% % aa = fscanf(dp821A);
% % display(aa);
% fprintf(dp821A,':OUTP:OCP:VAL CH2,1.500');
% fprintf(dp821A,':OUTP:OCP:VAL? CH2');
% aa = fscanf(dp821A);
% display(aa);
% fprintf(dp821A,':OUTP:OCP CH2,ON');
% fprintf(dp821A,':OUTP CH2,ON');
% fprintf(dp821A,':APPL CH1,24,1.000');
% fprintf(dp821A,':OUTP CH1,ON');

bias_name = 780;
pause(0.5);

M = 8;
data_length_initial = 10000;
zero_length = 3000;
zero_length_forsyn = 200;
ls_order = 50;
num_of_windows = 100;

times = 6;
origin_rate_tmp = 10e6;
f_rate = 160e6;
d_rate_tmp = origin_rate_tmp*times;

filter_order = 2000;     % filter order used in function sam_rate_con
rp = 0.00057565;      
rst = 1e-4;       % filter parameter used in function sam_rate_con

%% Generate filter
[origin_rate , d_rate , upf_transmit , dof_transmit , filter_transmit , upf_receive , dof_receive , filter_receive] = ruo_filter_gen(origin_rate_tmp , f_rate , times , filter_order , rp , rst);
rate_change = 100*abs(origin_rate-origin_rate_tmp)/origin_rate_tmp;
fprintf("rate's change = %.8f%% \n",rate_change);
ups_time = upf_transmit/dof_transmit;
% data_length = round(data_length_initial*8/ups_time);
data_length = data_length_initial;

filter_10m = firceqrip(filter_order,0.5e6/(origin_rate/2),[rp rst],'high');
%% Generate pilot
% pilot = pilot_gen([0 0 0 1 1 1 1]);
% pilot = pilot_gen([1 1 1 0 0 0 1 1 1]);
pilot_ini = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 0 1]);    %  2^11 = 2048
% pilot = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 1 0 1 1]);   %  2^13 = 8192
pilot_bpsk_ini = pilot_ini*2-1;
pilot_bpsk_tmp = conv(pilot_bpsk_ini,filter_10m);
pilot_bpsk_tmp = pilot_bpsk_tmp((length(filter_10m)+1)/2 : length(pilot_bpsk_tmp)-(length(filter_10m)-1)/2);
pilot_bpsk_tmp(pilot_bpsk_tmp <= 0) = -1;
pilot_bpsk_tmp(pilot_bpsk_tmp > 0) = 1;
pilot_bpsk = pilot_bpsk_tmp;
pilot = ruo_pamdemod(pilot_bpsk,2);
pilot_length = length(pilot);

%%
% origin_rate = gpuArray(double(origin_rate));
% f_rate = gpuArray(double(f_rate));
% d_rate = gpuArray(double(d_rate));
% num_of_windows = gpuArray(double(num_of_windows));
% times = gpuArray(double(times));
% pilot_length = gpuArray(double(pilot_length));
% delay = gpuArray(double(delay));
% filter_transmit = gpuArray(double(filter_transmit));
% filter_receive = gpuArray(double(filter_receive));
% upf_transmit = gpuArray(double(upf_transmit));
% dof_transmit = gpuArray(double(dof_transmit));
% upf_receive = gpuArray(double(upf_receive));
% dof_receive = gpuArray(double(dof_receive));
% h_channel = gpuArray(double(h_channel));
% h_channel_delay = gpuArray(double(h_channel_delay));

%%
t = datetime('now');
bias = 0.3;

% amp_begin = -4;
% amp_end = 60;
% amp_inf = 40;
amp_begin = 0.0015;
amp_end = 1;
amp_step = (amp_end-amp_begin)/50;
amp_inf = 0.6;
total_num = (amp_end-amp_begin)/amp_step+1;

fprintf('rand data \n');
amp = 1;

looptime = 0;
save_num = 50;

snr = 0;
while(looptime < save_num)
    
    looptime = looptime+1;
    %%
%     data_mpam_ini = [];
%     while size(data_mpam_ini) < data_length
%         x=rand*2-1;
%         y=rand;
%         if (y < x^2)
%             data_mpam_ini = [data_mpam_ini x];
%         end
%     end
    data_ini = randi([0,1],[1,data_length]);
    data_mpam_ini = real(pammod(data_ini,2));
        
    if looptime == 1
        data0 = data_mpam_ini;
    else
        data0_tmp = data_mpam_ini;
        data0 = [data0 data0_tmp];
    end
 
%     data_mpam_ini = rand(1,data_length);
    data_mpam_tmp = conv(data_mpam_ini,filter_10m);
    data_mpam = data_mpam_tmp((length(filter_10m)+1)/2 : length(data_mpam_tmp)-(length(filter_10m)-1)/2);
    data_mpam = data_mpam./norm(data_mpam,2)*sqrt(length(data_mpam));
    data = data_mpam;
    if looptime == 1
       data1 = data;
    else
       data1_tmp = data;
       data1 = [data1 data1_tmp];
    end
    %%
    signal_upsample = upsample(data,upf_transmit);  % upsample signal
    signal_upsample = signal_upsample*upf_transmit;
    
    signal_passfilter = conv(signal_upsample,filter_transmit);                        % Filter the upsampled signal
    signal_passfilter = signal_passfilter((length(filter_transmit)+1)/2:length(signal_passfilter)-(length(filter_transmit)-1)/2);
    if looptime == 1
        data2 = signal_passfilter;
    else
        data2_tmp = signal_passfilter;
        data2 = [data2 data2_tmp];
    end
    
    if fix(length(signal_passfilter)/2)+500*dof_transmit <= length(signal_passfilter)
        signal_fil_judge = signal_passfilter(1,fix(length(signal_passfilter)/2):fix(length(signal_passfilter)/2)+500*dof_transmit);
    else
        signal_fil_judge = signal_passfilter(1,fix(length(signal_passfilter)/2):end);
    end
    energy = zeros(1,dof_transmit);
    for downsample_phase = 0:dof_transmit-1
        sampled_signal_judge = downsample(signal_fil_judge,dof_transmit,downsample_phase);
        energy(downsample_phase+1) = norm(sampled_signal_judge,2);
    end
    [~,judge_phase] = max(energy);
    temporary_phase = fix(length(signal_passfilter)/2)+judge_phase-1;
    real_phase = mod(temporary_phase,dof_transmit);
    if real_phase == 0
        real_phase = dof_transmit;
    end                                                  % Find the downsample phase with max energy
    
    signal_final_sampled = downsample(signal_passfilter,dof_transmit,real_phase-1);
    if looptime == 1
        data3 = signal_final_sampled;
    else
        data3_tmp = signal_final_sampled;
        data3 = [data3 data3_tmp];
    end
end
hist(data0,100);
figure
hist(data1,100);
figure
hist(data2,100);
% figure
% hist(data3,100);