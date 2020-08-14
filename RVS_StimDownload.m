function dur  = RVS_stimDownload(stim_mat,pin_map,dt)
% This function downloads the stimulus to the QuaeroSys 
% Independent of the QuaeroBox
% after completion of the download the time of download is handed back

start_time = time;

% possible values of a single pin range from 0 to 4095
maxhub = 4095;
nsamples = size(stim_mat,2);

calllib('stimlib0', 'setPinBlock', 0, 5, pin_map(1,2), pin_map(1,1), pin_map(2,2), pin_map(2,1), pin_map(3,2), pin_map(3,1), pin_map(4,2), pin_map(4,1));
calllib('stimlib0', 'setPinBlock', 1, 5, pin_map(1,4), pin_map(1,3), pin_map(2,4), pin_map(2,3), pin_map(3,4), pin_map(3,3), pin_map(4,4), pin_map(4,3));


for n = 1:(nsamples-1)
    calllib('stimlib0', 'setDAC', 1, 0);
    calllib('stimlib0', 'setDAC', 2, round(stim_mat(n)*maxhub));    
    calllib('stimlib0', 'wait', 5, dt/0.5); %
end
 
dur = time-start_time;

