% BMED 7610 HW 2 

clear
close all hidden
% read matrix file
cortical_spikes = load('HWK2_prob1A_data.mat','-mat','spike_times');

t = 1:9000;
dt = 100; %bin size for the window in msec

% a) plot the spike train of the first trial
spk_times = cortical_spikes.spike_times;
% figure(1)
% plot(spk_times(1, :), 1:36)
% title('The spike train of the first trial')
% xlabel('time (msec)')
% ylabel('spike number');
% b) a raster plot of all the trials

x = spk_times(1, :);
y = ones([1 36]);
% plot_raster(x, y);
% figure(2)
hold on
for trial = 1:28
    x = spk_times(trial, :);
    y = trial*ones([1 36]);
%     plot_raster(x, y);
end
hold off
% title('The spike train of trials')
% xlabel('time (msec)')
% ylabel('trials');
%% c) calculate firing rate of the first trial using fixed bins of 100 msec

firing_rate = zeros([1 8900]);
num_spikes = zeros([1 8900]);
for time = 100:9000
    num_spikes(time) = sum(spk_times(1, :)>time & spk_times(1, :)<(time+dt));
    firing_rate(time) = num_spikes(time)/dt;
    if mod(time, 100) == 0
        num_spikes(time-99:time)= num_spikes(time);
        firing_rate(time-99:time) = firing_rate(time);
        time = time + dt;
    end
end
firing_rate = firing_rate*10; %converting to Hz
%% d) Calculate firing rate of the first trial using a sliding rectangular bin of 100 msec

firing_rate_sliding = zeros([1 9000]);
num_spikes_sliding = zeros([1 9000]);
for time=1:8900
    num_spikes_sliding(time)=sum(spk_times(1, :)>time & spk_times(1, :)<(time+dt));
    firing_rate_sliding(time) = num_spikes_sliding(time)/dt;
end
firing_rate_sliding = firing_rate_sliding*10; %converting to Hz
%% e) Calculate firing rate of the first trial by convolution of a gaussian window of width sigma = 100 msec

firing_rate_gaus = zeros([1 9000]);
for time = 1:9000
    firing_rate_gaus(time) = 0;
    for j=1:36
        firing_rate_gaus(time) = firing_rate_gaus(time) + w_g(time - spk_times(1, j), dt);
    end
end
firing_rate_gaus = firing_rate_gaus*10; %converting to Hz
%% f)repeat part a-e with dt of 10 msec

dt = 10;
firing_rate_ten = zeros([1 8900]);
num_spikes_ten = zeros([1 8900]);
for time = 100:9000
    num_spikes_ten(time) = sum(spk_times(1, :)>time & spk_times(1, :)<(time+dt));
    firing_rate_ten(time) = num_spikes_ten(time)/dt;
    if mod(time, 10) == 0
        num_spikes_ten(time-9:time)= num_spikes_ten(time);
        firing_rate_ten(time-9:time) = firing_rate_ten(time);
        time = time + dt; 
    end
end
firing_rate_ten = firing_rate_ten*10; %converting to Hz
firing_rate_sliding_ten = zeros([1 9000]);
num_spikes_sliding_ten = zeros([1 9000]);
for time=1:8900
    num_spikes_sliding_ten(time)=sum(spk_times(1, :)>time & spk_times(1, :)<(time+dt));
    firing_rate_sliding_ten(time) = num_spikes_sliding_ten(time)/dt;
end
firing_rate_sliding_ten = firing_rate_sliding_ten*10; %converting to Hz
firing_rate_gaus_ten = zeros([1 9000]);
for time = 1:9000
    firing_rate_gaus_ten(time) = 0;
    for j=1:36
        firing_rate_gaus_ten(time) = firing_rate_gaus_ten(time) + w_g(time - spk_times(1, j), dt);
    end
end
firing_rate_gaus_ten = firing_rate_gaus_ten*10; %converting to Hz
%% Plots
figure(1)
subplot(3,1,1)
plot(t, firing_rate);
title('The first trial firing rate using fixed bins of 100 msec')
xlabel('time (msec)')
ylabel('rate (Hz)')

subplot(3,1,2)
plot(t, firing_rate_sliding)
title('The firing rate of first trail using sliding bins of width 100 msec.')
xlabel('time (msec)')
ylabel('rate (Hz)')

subplot(3,1,3)
plot(t, firing_rate_gaus)
title('The firing rate of first trail using gaussian with sigma of 100 msec.')
xlabel('time (msec)')
ylabel('rate (Hz)')

figure(2)
subplot(3,1,1)
plot(t, firing_rate_ten);
title('The first trial firing rate using fixed bins of 10 msec')
xlabel('time (msec)')
ylabel('rate (Hz)')

subplot(3,1,2)
plot(t, firing_rate_sliding_ten)
title('The firing rate of first trail using sliding bins of width 10 msec.')
xlabel('time (msec)')
ylabel('rate (Hz)')

subplot(3,1,3)
plot(t, firing_rate_gaus_ten)
title('The firing rate of first trail using gaussian with sigma of 10 msec.')
xlabel('time (msec)')
ylabel('rate (Hz)')

%% g) the firing rate for the average of all trials, using bins of 100 and 
%     10 msec.

dt = 100;
firing_rate_avg = zeros([1 8900]);
num_spikes_avg = zeros([1 8900]);
for time = 100:9000
    num_spikes_avg(time) = mean(sum(spk_times(:, :)>time & spk_times(:, :)<(time+dt))/28);
    firing_rate_avg(time) = num_spikes_avg(time)/dt;
    if mod(time, 100) == 0
        num_spikes_avg(time-99:time)= num_spikes_avg(time);
        firing_rate_avg(time-99:time) = firing_rate_avg(time);
        time = time + dt;
    end
end
firing_rate_avg = firing_rate_avg*10; %converting to Hz

dt = 10;
firing_rate_avg_ten = zeros([1 8900]);
num_spikes_avg_ten = zeros([1 8900]);
for time = 100:9000
    num_spikes_avg_ten(time) = mean(sum(spk_times(:, :)>time & spk_times(:, :)<(time+dt))/28);
    firing_rate_avg_ten(time) = num_spikes_avg_ten(time)/dt;
    if mod(time, 10) == 0
        num_spikes_avg_ten(time-9:time)= num_spikes_avg_ten(time);
        firing_rate_avg_ten(time-9:time) = firing_rate_avg_ten(time);
        time = time + dt;
    end
end
firing_rate_avg_ten = firing_rate_avg_ten*10; %converting to Hz

figure(3)
subplot(2,1,1)
plot(t, firing_rate_avg);
title('The avg firing rate for all trials using fixed bins of 100 msec')
xlabel('time (msec)')
ylabel('rate (Hz)')

subplot(2,1,2)
plot(t, firing_rate_avg_ten);
title('The avg firing rate for all trials using fixed bins of 10 msec')
xlabel('time (msec)')
ylabel('rate (Hz)')
%% h) the firing rate, calculated by convolving with a rectangular window 
%      over each trial separately, and subsequently averaging across trials

win100 = ones(100,1);
win10 = ones(10,1);
for trial=1:28
    firing_rate_conv(trial, :) = conv(spk_times(trial, :), win100);
    firing_rate_conv_ten(trial, :) = conv(spk_times(trial, :), win10);
    
end
firing_rate_conv = firing_rate_conv*10; %converting to Hz
firing_rate_conv_ten = firing_rate_conv_ten*10; %converting to Hz
avg_firing_rate_conv = mean(firing_rate_conv);
avg_firing_rate_conv_ten = mean(firing_rate_conv_ten);

%% Plot for part h

figure(3)
subplot(2,1,1)
plot(avg_firing_rate_conv);
title('The avg firing rate for all trials using conv of 100 msec window')
xlabel('time (msec)')
ylabel('rate (Hz)')

subplot(2,1,2)
plot(avg_firing_rate_conv_ten);
title('The avg firing rate for all trials using conv of 10 msec window')
xlabel('time (msec)')
ylabel('rate (Hz)')

%% Gausian Window Function

function W_g = w_g(t,sig_w)
    W_g = (1/(sqrt(2*pi)*sig_w))*exp(-t.^2/(2*sig_w.^2));
end
%% Window function
