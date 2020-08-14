function sequences = RVSsequences3(Stimuli,sequence_length)


fre_sequences = [1 3 2;
                 2 4 3;
                 1 4 3;
                 2 4 1];
fre_rand = [];
for i = 1:size(fre_sequences)
    fre_rand = [perms(fre_sequences(i,:));fre_rand];  
end        
        
        

sequences.fre_sequences = fre_rand;