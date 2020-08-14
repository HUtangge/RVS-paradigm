function Stimuli = RVS_Stimuli3()

% Weber-Fechner Korrektur

n_stim = 8;
lowest = 12;
dist1  = 2.5;
const = dist1/lowest;
freqs = [lowest];
for i = 1:(n_stim-1)
    next = freqs(i)*const+freqs(i);
    freqs = [freqs next];
end

%plot(1:n_stim,freqs);
Stimuli.Frequencies = freqs;
Stimuli.patfreq = mean(Stimuli.Frequencies);

Stimuli.FreqTarLow = Stimuli.Frequencies(1:3);

Stimuli.FreqTarHigh = Stimuli.Frequencies(6:8);

Stimuli.Pattern{1}= [ 1 1 1 1;
                      2 2 1 1;
                      2 2 1 1;
                      2 2 1 1];

Stimuli.Pattern{2}= [ 1 1 1 1 ;
                      1 1 1 1;
                      2 2 2 1;
                      2 2 2 1];

Stimuli.Pattern{3}= [ 2 2 1 1
                      2 2 1 1;
                      2 2 1 1;
                      1 1 1 1];
  
Stimuli.Pattern{4}= [ 2 2 2 1;
                      2 2 2 1;
                      1 1 1 1;
                      1 1 1 1];

Stimuli.Pattern{5}= [ 1 1 2 2;
                      1 1 2 2;
                      1 1 2 2;
                      1 1 1 1];

Stimuli.Pattern{6}= [ 1 2 2 2;
                      1 2 2 2;
                      1 1 1 1;
                      1 1 1 1];  

Stimuli.Pattern{7}= [ 1 1 1 1;
                      1 1 2 2;
                      1 1 2 2;
                      1 1 2 2];

Stimuli.Pattern{8}= [ 1 1 1 1;
                      1 1 1 1;
                      1 2 2 2;
                      1 2 2 2];  



Stimuli.PatTarHigh{1} = Stimuli.Pattern{1};
Stimuli.PatTarHigh{2} = Stimuli.Pattern{2};
Stimuli.PatTarHigh{3} = Stimuli.Pattern{7};
Stimuli.PatTarHigh{4} = Stimuli.Pattern{8};

Stimuli.PatTarLow{1} = Stimuli.Pattern{3};
Stimuli.PatTarLow{2} = Stimuli.Pattern{4};
Stimuli.PatTarLow{3} = Stimuli.Pattern{5};
Stimuli.PatTarLow{4} = Stimuli.Pattern{6};




                  
                  
                  