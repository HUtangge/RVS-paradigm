fre_sequences = [perms([3,2,1]);
    perms([4,2,1]);
    perms([4,3,1]);
    perms([4,3,2])];

%% Randomize the rehearsal sequence
trial_per_block = randperm(length(fre_sequences));

%% Assign them rehearsal duration
% group = find(fre_sequences == 1);
group = [];

% subgroup means four groups with all permutation
subgroup = [];
for j = 1:4
    group = [group find(fre_sequences == j)];
    for i = 1:length(group)
        if group(i,j) > 24 & group(i,j) <= 48
            group(i,j) = group(i,j) - 24;
        elseif group(i,j) > 48
            group(i,j) = group(i,j) - 48;
        else
            group(i,j) = group(i,j);
        end
    end
    subgroup(:,j) = trial_per_block(find(ismember(trial_per_block, group(:,j)) == 0)).';
    subgroup_size = size(subgroup);
    
    % assign rehearsal duration 8s to two sequences
    subgroup(randperm(subgroup_size(1),2),j+4+(j-1)*3) = 8;
    subgroup = [subgroup fre_sequences(subgroup(:,j),:)];
    dur_8 = [-1 -1];
    rand_10 = 0;
    dur_8 = subgroup(find(subgroup(:,j+4+(j-1)*3) == 8),j+7+(j-1)*3);
    dur_10_idx = find(ismember(subgroup(:,j+5+(j-1)*3),dur_8));
    dur_10_idx = setdiff(dur_10_idx, find(subgroup(:,j+4+(j-1)*3) == 8));
    while 1
        %ismember(subgroup(:,6),subgroup(find(subgroup(:,5) == 8),8));
        % randomly assign 10s to the rest 4 sequences and make sure they
        % are the same as the last frequency assigned to 8s rehearsal
        % duration
        rand_10_idx = dur_10_idx(randperm(length(dur_10_idx),2));
        rand_10 = subgroup(rand_10_idx, j+5+(j-1)*3);
        if sum(dur_8 == rand_10) == 2
            break
        elseif sum(dur_8 == rand_10) == 0
            break
        end           
    end
    subgroup(rand_10_idx,j+4+(j-1)*3) = 10;
    
    % the rest two frequences have a 6s rehearsal duration
    subgroup(find(subgroup(:,j+4+(j-1)*3) == 0), j+4+(j-1)*3) = 6;
    
end

%% Assign stimuli duration and stimuli sequence
sub_group = {};
for i = 1:4
    sub_group = {};
    test = [(subgroup(:,i+5+(i-1)*3) - subgroup(:,i+6+(i-1)*3)) (subgroup(:,i+6+(i-1)*3) - subgroup(:,i+7+(i-1)*3))];    
    for m = 1:2
        sub_group{end+1} = strcat('subgroup_', num2str(m));
    end   
    sub_group{1} = find((test(:,1) < 0 & test(:,2) < 0) | test(:,1) == max(max(test))| test(:,2) == max(max(test)));
    sub_group{2} = setdiff(find(subgroup(:,i+5+(i-1)*3) ~= 0), sub_group{1});
    for n = 1:2
        stim_dur = [4 6 8];
        first_letter_1 = 10;
        for j = 1:length(sub_group{n})
            part_stim_dur = stim_dur(randperm(length(stim_dur),1));
            subgroup(sub_group{n}(j),20+i+(i-1)*3) = part_stim_dur;          
            switch subgroup(sub_group{n}(j),20+i+(i-1)*3)
                case 4
                    first_letter = subgroup(sub_group{n}(j),i+7+(i-1)*3);
                    second_letter = subgroup(sub_group{n}(j),i+5+(i-1)*3);
                    third_letter = subgroup(sub_group{n}(j),i+6+(i-1)*3);
                    %stim_seq = [first_letter subgroup(sub_group{n}(j),i+5+(i-1)*3) subgroup(sub_group{n}(j),i+7+(i-1)*3)];
                case 6
                    first_letter = subgroup(sub_group{n}(j),i+5+(i-1)*3);
                    second_letter = subgroup(sub_group{n}(j),i+6+(i-1)*3);
                    third_letter = subgroup(sub_group{n}(j),i+7+(i-1)*3);
                    %stim_seq = [first_letter subgroup(sub_group{n}(j),i+6+(i-1)*3) subgroup(sub_group{n}(j),i+5+(i-1)*3)];
                case 8
                    first_letter = subgroup(sub_group{n}(j),i+6+(i-1)*3);
                    second_letter = subgroup(sub_group{n}(j),i+7+(i-1)*3);
                    third_letter = subgroup(sub_group{n}(j),i+5+(i-1)*3);
                    %stim_seq = [first_letter subgroup(sub_group{n}(j),i+7+(i-1)*3) subgroup(sub_group{n}(j),i+6+(i-1)*3)];
            end

            if first_letter_1 == first_letter;
                subgroup(sub_group{n}(j),20+i+(i-1)*3) = setdiff(stim_dur, part_stim_dur);
                first_letter = 0;
                second_letter = 0;
                third_letter = 0;
            end
            subgroup(sub_group{n}(j),20+i+(i-1)*3+1) = first_letter;
            subgroup(sub_group{n}(j),20+i+(i-1)*3+2) = second_letter;
            subgroup(sub_group{n}(j),20+i+(i-1)*3+3) = third_letter;
            first_letter_1 = first_letter;            
            
            stim_dur(stim_dur == subgroup(sub_group{n}(j),20+i+(i-1)*3)) = [];
            
        end
    subgroup(sub_group{n}(find(subgroup(sub_group{n},20+i+(i-1)*3+1) == 0)), 20+i+(i-1)*3+1) = setdiff(subgroup(sub_group{n}(j),(i+5+(i-1)*3): (i+7+(i-1)*3)), subgroup(sub_group{n},20+i+(i-1)*3+1).'); 
    subgroup(sub_group{n}(find(subgroup(sub_group{n},20+i+(i-1)*3+2) == 0)), 20+i+(i-1)*3+2) = setdiff(subgroup(sub_group{n}(j),(i+5+(i-1)*3): (i+7+(i-1)*3)), subgroup(sub_group{n},20+i+(i-1)*3+2).'); 
    subgroup(sub_group{n}(find(subgroup(sub_group{n},20+i+(i-1)*3+3) == 0)), 20+i+(i-1)*3+3) = setdiff(subgroup(sub_group{n}(j),(i+5+(i-1)*3): (i+7+(i-1)*3)), subgroup(sub_group{n},20+i+(i-1)*3+3).'); 
    
    end
end


%% Combine all
% column 1 2 3 is stimuli sequence, column 4 is stimuli duration, 
% column 5 6 7 is rehearsal sequence, column 8 is rehearsal duration
for i = 1:4
    fre_sequences(subgroup(:,i), 4) = subgroup(:,i+4+(i-1)*3);
    fre_sequences(subgroup(:,i), 5:7) = subgroup(:,(20+i+(i-1)*3+1):(20+i+(i-1)*3+3));
    fre_sequences(subgroup(:,i), 8) = subgroup(:,20+i+(i-1)*3);    
end
fre_sequences = fre_sequences(:,[5 6 7 8 1 2 3 4]);

%% Some statistical stuff
%a = [fre_sequences(find(fre_sequences(:,8) == 8),5) fre_sequences(find(fre_sequences(:,8) == 8),6) fre_sequences(find(fre_sequences(:,8) == 10),5)];
%b = [fre_sequences(find(fre_sequences(:,8) == 8),5), fre_sequences(find(fre_sequences(:,8) == 10),5)];
%c = fre_sequences(:,5);
%d = fre_sequences(find(fre_sequences(:,8) == 8),5:7)
e = find(fre_sequences(:,4) == 4 & fre_sequences(:,8) == 8)
%length(find(b == 1))

%rand_fre(:,:,1) = fre_sequences(:,[1 2 3 4]);
%rand_fre(:,:,2) = fre_sequences(:,[5 6 7 8]);
%permute(rand_fre, [3 1 2]);





