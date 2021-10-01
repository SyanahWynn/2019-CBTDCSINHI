%% OPEN .txt FILE
clear all 
close all
% turn off unicode warning
warning off

% set researcher
prompt          = '\nWho are you?\n';
researcher      = input(prompt,'s');

clear prompt

% variables
if strcmp('researcher',researcher) || strcmp('Researcher',researcher)
    inputloc  = '???';
    outputloc = '???'; 
end
outputfile = fullfile(outputloc,'RawData');
stims = [0,1;0,1;0,1;0,1;1,0;0,1;1,0;1,0;1,0;1,0;0,1;1,0;0,1;1,0;1,0;1,0;0,1;0,1;0,1;1,0;1,0;1,0;0,1;0,1;1,0;0,1];


if ~exist(outputloc, 'dir')
  mkdir(outputloc);
end

% location of .txt file
txtdir      = fullfile(inputloc,'*txt*');
txtdf       = dir(txtdir);
txtfiles    = {txtdf.name};
clear txtdir txtdf

preproc = true; % to do or to not to 'preprocessing'
pre_ana = true; % to do or to not to the summarizing of the data

% to loop or not to loop
prompt          = '\nWho do you want to analyse? (all, 101, 102, etc)\n';
strt      = input(prompt,'s');
clear prompt
% if all participants have to be analysed, create a loop
% start loop or determine current ppn
if strcmp(strt,'all')
    fprintf('Loop true')
    strt=0; % make a not char thing
else
    fprintf('Loop false')
end

%%%%%%%%%%%%%%%%%%%
%% PREPROCESSING %%
%%%%%%%%%%%%%%%%%%%
if preproc

    %%%%%%%%%%%%%%%%%%%%
    %% LOAD THE DATA %%%
    %%%%%%%%%%%%%%%%%%%%

    for t=1:length(txtfiles)
        % get the task name 
        ind = find(ismember(txtfiles{t},'_'));
        cur_task = txtfiles{t}(1:ind-1);
        % get the participant number
        ind = find(ismember(txtfiles{t},'-'));
        cur_ppn = txtfiles{t}(ind(1)+1:ind(2)-1);
        if cur_ppn == strt || ~ischar(strt)
            % get the session number
            cur_ses = txtfiles{t}(ind(2)+1);
            fprintf('########\nLoading... %s task of participant %s during session %s\n########',cur_task,cur_ppn,cur_ses)
            clear ind*
            % get the file identifier
            % add the unicode to correct for 16 bit/8 bit confusion
            fid = fopen(fullfile(inputloc,txtfiles{t}), 'r','n','Unicode');
            % get the first line of the file, for the while loop
            line = fgetl(fid);
            % open an empty array to hold the data from the file
            full_txt = {};
            % while there is data is the file, keep going
            while ischar(line)
                % add the line to the array
                full_txt{end+1} = line;
                % get the next line from the file
                line = fgetl(fid);
            end
            fclose(fid);
            clear line fid txt_dir

            %% GET HEADER VARIABLES

            % loop over the array
            hdr_start = false;
            hdr_end = false;
            for i=1:length(full_txt)
                % find the start of the header information
                if all(ismember('Header Start',full_txt{i}))
                    hdr_start = i;
                end
                % find the end of the header information
                if all(ismember('Header End',full_txt{i}))
                    hdr_end = i;
                end
                % if both are found exit loop
                if hdr_start && hdr_end
                    break
                end
            end
            % create table to store data
            Data_temp = table;
            % store the header variables in the table
            for i=hdr_start+1:hdr_end-1
                % find the colon
                colon = find(ismember(full_txt{i},':'));
                var_name = full_txt{i}(1:colon-1);
                value = full_txt{i}(colon+2:end);
                % remove dots from variable name
                var_name=var_name(~ismember(var_name,'.'));
                % check for empty values
                if isempty(value)
                    value = 'None';
                end
                % check if variable is a number
                if all(ismember(value, '0123456789.eEdD'))
                    % add variable to table (as number)
                    value = str2double(value);
                    evalc(sprintf('Data_temp.%s=%f;',var_name,value));
                else
                    % add variable to table (as string)
                    evalc(sprintf('Data_temp.%s=''%s'';',var_name,value));
                end
            end
            clear hdr*

            %% GET THE OTHER VARIABLES

            % loop over the array
            logfr_start = [];
            logfr_end = [];
            for i=1:length(full_txt)
                % find the starts of the LogFrame (start of trial)
                if all(ismember('LogFrame Start',full_txt{i}))
                    logfr_start = [logfr_start i];
                end
                % find the ends of the LogFrame (end of trial)
                if all(ismember('LogFrame End',full_txt{i}))
                    logfr_end = [logfr_end i];
                end
            end
            % check whether the amount of starts matches the amount of ends
            if isequal(length(logfr_start),length(logfr_end))
                % loop over the trials
                for i=1:length(logfr_start)
                    % store the variables in the table
                    Data_temp.trial(i,1)=i;
                    for r=logfr_start(i)+1:logfr_end(i)-1
                        % find the colon
                        colon = find(ismember(full_txt{r},':'));
                        var_name = full_txt{r}(1:colon-1);
                        value = full_txt{r}(colon+2:end);
                        % remove dots from variable name
                        var_name=var_name(~ismember(var_name,'.'));
                        % check for empty values
                        if isempty(value)
                            value = 'None';
                        end
                        % check if variable is a number
                        if all(ismember(value, '0123456789.eEdD'))
                            % add variable to table (as number)
                            value = str2double(value);
                            evalc(sprintf('Data_temp.%s(i,1)=%f;',var_name,value));
                        else
                            % add variable to table (as string)
                            evalc(sprintf('Data_temp.%s(i,1:%d)=''%s'';',var_name,length(value),value));
                        end
                        clear value var_name colon
                    end
                end
            else
                error('the amount of LogFrame starts does not match the amount of ends. Check if the data file is not intact')
            end
            clear logfr* r i
            evalc(sprintf('Data.%s.pp_%s.session_%s = Data_temp;',cur_task,cur_ppn,cur_ses));
            clear Data_temp full_txt

        end
    end
    clear t

    % save the data
    save(outputfile,'Data')
    outputfile2 = strcat(outputfile,'_',date);
    save(outputfile2,'Data')
    fprintf('DATA SAVED')
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SUMMARY DATA PER PARTICIPANT (AND PER SESSION) %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if pre_ana
    
    % load the data
    load(outputfile)
    
    %% Temporal Discounting Task

    % VARIABLES
    delayed_reward = 10;
    delays = [0 2,14,30,180,365]; % add 0 at the start to calculate the are of the first trapeziod
    delay_norm = delays./delays(end);

    % loop over participants
    ppns = fieldnames(Data.TD); 
    ppns=ppns(contains(ppns,'pp'));
    % if all participants have to be analysed, create a loop
    % start loop or determine current ppn
    if ~ischar(strt)
        fprintf('Loop true')
        strt = 1:length(ppns);
    else
        fprintf('Loop false')
        strt=1;
    end
    for p=strt
        % loop over the sessions
        sess = eval(sprintf('fieldnames(Data.TD.%s)',ppns{p}));
        for s=1:length(sess)
            % store the currently relevant data in a temp variable
            cur_data = eval(sprintf('Data.TD.%s.%s',ppns{p},sess{s}));
            % get the variables of interest for this dataset
            ppn.Subject    = cur_data.Subject(1);
            ppn.Session    = cur_data.Session(1);
            ppn.Age        = cur_data.Age(1);
            ppn.Sex        = cur_data.Sex(1);
            ppn.SV_2       = cur_data.SV(logical(sum(ismember(cur_data.DRdelay,'2'),2)==1));
            ppn.SV_14      = cur_data.SV(logical(sum(ismember(cur_data.DRdelay,'14'),2)==2));
            ppn.SV_30      = cur_data.SV(logical(sum(ismember(cur_data.DRdelay,'30'),2)==2));
            ppn.SV_180     = cur_data.SV(logical(sum(ismember(cur_data.DRdelay,'180'),2)==3));
            ppn.SV_365     = cur_data.SV(logical(sum(ismember(cur_data.DRdelay,'365'),2)==3));
            ppn.RT_Choice1 = cur_data.Choice1RT(cur_data.Choice1RT~=0);
            ppn.RT_Choice2 = cur_data.Choice2RT(cur_data.Choice2RT~=0);
            ppn.RT_Choice3 = cur_data.Choice3RT(cur_data.Choice3RT~=0);
            ppn.RT_Choice4 = cur_data.Choice4RT(cur_data.Choice4RT~=0);
            ppn.RT_Choice5 = cur_data.Choice5RT(cur_data.Choice5RT~=0);
            ppn.RT_Choice6 = cur_data.Choice6RT(cur_data.Choice6RT~=0);
            ppn.RT_2       = [cur_data.Choice1RT(logical(sum(ismember(cur_data.DRdelay,'2'),2)==1));
                                cur_data.Choice2RT(logical(sum(ismember(cur_data.DRdelay,'2'),2)==1));
                                cur_data.Choice3RT(logical(sum(ismember(cur_data.DRdelay,'2'),2)==1));
                                cur_data.Choice4RT(logical(sum(ismember(cur_data.DRdelay,'2'),2)==1));
                                cur_data.Choice5RT(logical(sum(ismember(cur_data.DRdelay,'2'),2)==1));
                                cur_data.Choice6RT(logical(sum(ismember(cur_data.DRdelay,'2'),2)==1))];
            ppn.RT_14       = [cur_data.Choice1RT(logical(sum(ismember(cur_data.DRdelay,'14'),2)==2));
                                cur_data.Choice2RT(logical(sum(ismember(cur_data.DRdelay,'14'),2)==2));
                                cur_data.Choice3RT(logical(sum(ismember(cur_data.DRdelay,'14'),2)==2));
                                cur_data.Choice4RT(logical(sum(ismember(cur_data.DRdelay,'14'),2)==2));
                                cur_data.Choice5RT(logical(sum(ismember(cur_data.DRdelay,'14'),2)==2));
                                cur_data.Choice6RT(logical(sum(ismember(cur_data.DRdelay,'14'),2)==2))];
            ppn.RT_30       = [cur_data.Choice1RT(logical(sum(ismember(cur_data.DRdelay,'30'),2)==2));
                                cur_data.Choice2RT(logical(sum(ismember(cur_data.DRdelay,'30'),2)==2));
                                cur_data.Choice3RT(logical(sum(ismember(cur_data.DRdelay,'30'),2)==2));
                                cur_data.Choice4RT(logical(sum(ismember(cur_data.DRdelay,'30'),2)==2));
                                cur_data.Choice5RT(logical(sum(ismember(cur_data.DRdelay,'30'),2)==2));
                                cur_data.Choice6RT(logical(sum(ismember(cur_data.DRdelay,'30'),2)==2))];               
            ppn.RT_180       = [cur_data.Choice1RT(logical(sum(ismember(cur_data.DRdelay,'180'),2)==3));
                                cur_data.Choice2RT(logical(sum(ismember(cur_data.DRdelay,'180'),2)==3));
                                cur_data.Choice3RT(logical(sum(ismember(cur_data.DRdelay,'180'),2)==3));
                                cur_data.Choice4RT(logical(sum(ismember(cur_data.DRdelay,'180'),2)==3));
                                cur_data.Choice5RT(logical(sum(ismember(cur_data.DRdelay,'180'),2)==3));
                                cur_data.Choice6RT(logical(sum(ismember(cur_data.DRdelay,'180'),2)==3))];
            ppn.RT_365       = [cur_data.Choice1RT(logical(sum(ismember(cur_data.DRdelay,'365'),2)==3));
                                cur_data.Choice2RT(logical(sum(ismember(cur_data.DRdelay,'365'),2)==3));
                                cur_data.Choice3RT(logical(sum(ismember(cur_data.DRdelay,'365'),2)==3));
                                cur_data.Choice4RT(logical(sum(ismember(cur_data.DRdelay,'365'),2)==3));
                                cur_data.Choice5RT(logical(sum(ismember(cur_data.DRdelay,'365'),2)==3));
                                cur_data.Choice6RT(logical(sum(ismember(cur_data.DRdelay,'365'),2)==3))];
            % RTs
            ppn.RT_Mean    = mean([ppn.RT_Choice1;ppn.RT_Choice2;ppn.RT_Choice3;ppn.RT_Choice4;ppn.RT_Choice5;ppn.RT_Choice6]);
            ppn.RT_Choice1 = mean(ppn.RT_Choice1);
            ppn.RT_Choice2 = mean(ppn.RT_Choice2);
            ppn.RT_Choice3 = mean(ppn.RT_Choice3);
            ppn.RT_Choice4 = mean(ppn.RT_Choice4);
            ppn.RT_Choice5 = mean(ppn.RT_Choice5);
            ppn.RT_Choice6 = mean(ppn.RT_Choice6);
            ppn.RT_2       = mean(ppn.RT_2);
            ppn.RT_14      = mean(ppn.RT_14);
            ppn.RT_30      = mean(ppn.RT_30);
            ppn.RT_180     = mean(ppn.RT_180);
            ppn.RT_365     = mean(ppn.RT_365);

            % select the bonus
            % get the choices (now or delayed)
            now_later   = cur_data(1:end-1,(logical(strcmp('Choice1RESP',fieldnames(cur_data)) | ...
                                strcmp('Choice2RESP',fieldnames(cur_data)) | ...
                                strcmp('Choice3RESP',fieldnames(cur_data)) | ...
                                strcmp('Choice4RESP',fieldnames(cur_data)) | ...
                                strcmp('Choice5RESP',fieldnames(cur_data)) | ...
                                strcmp('Choice6RESP',fieldnames(cur_data)))));
            values      = cur_data(1:end-1,(logical(strcmp('AdAm1',fieldnames(cur_data)) | ...
                                strcmp('AdAm2',fieldnames(cur_data)) | ...
                                strcmp('AdAm3',fieldnames(cur_data)) | ...
                                strcmp('AdAm4',fieldnames(cur_data)) | ...
                                strcmp('AdAm5',fieldnames(cur_data)) | ...
                                strcmp('AdAm6',fieldnames(cur_data)))));
            SVs         = cur_data(1:end-1,(logical(strcmp('SV',fieldnames(cur_data)))) ...
                            );
            meaning_z   = {'IR','DR','DR','IR','IR','DR'}; % presented left on the screen. If participant choose this, 'z' was pressed
            % if pp is even, calculate bonus from session 2, of uneven calculate from session 1
            if ppn.Session==mod(ppn.Subject+1,2)+1
                % determine the delay 
                % exclude the 365 days option
                delays         = find(~logical(sum(ismember(cur_data(1:end-1,:).DRdelay,'365'),2)==3));
                % always choose the 3rd delay (delay is randomized in task)
                cur_delay = 3;
                % choose the choice based upon participant number
                pp_4 = [104,105,107,110,114,116,120,123,124];
                pp_6 = [101,102,103,111,112,117,118,121];
                pp_5 = [106,108,109,113,115,119,122,125,126];
                if ismember(ppn.Subject,pp_4) 
                    cur_choice = 4;
                elseif ismember(ppn.Subject,pp_6) 
                    cur_choice = 6;
                elseif ismember(ppn.Subject,pp_5) 
                    cur_choice = 5;
                end
                % check what the participant choose
                if strcmp(cell2mat(table2cell(now_later(delays(cur_delay),cur_choice))),'z') & strcmp(meaning_z(cur_choice),'IR') || ...
                    strcmp(cell2mat(table2cell(now_later(delays(cur_delay),cur_choice))),'m') & strcmp(meaning_z(cur_choice),'DR') 
                    % choice: IR
                    bonus_delay = cur_data(1,:).IRdelay;
                    bonus{1}      = sprintf('%.2f %s',cell2mat(table2cell(values(delays(cur_delay),cur_choice))),bonus_delay);
                elseif strcmp(cell2mat(table2cell(now_later(delays(cur_delay),cur_choice))),'m') & strcmp(meaning_z(cur_choice),'IR') || ...
                    strcmp(cell2mat(table2cell(now_later(delays(cur_delay),cur_choice))),'z') & strcmp(meaning_z(cur_choice),'DR') 
                    % choice: DR
                    bonus_delay    = cur_data(delays(cur_delay),:).DRdelay;
                    bonus{1}      = sprintf('%.2f %s',cell2mat(table2cell(SVs(delays(cur_delay),:))),bonus_delay);
                end
                %display(sprintf('\nbonus PP: %d\n%s',ppn.Subject, bonus{1}))
                ppn.bonus=bonus;
                clear pp_*
            else
                ppn.bonus{1}      = '';
            end

            %AUC (Myerson (2001) - Area under the curve as a measure of discounting)
            SV_norm = [delayed_reward ppn.SV_2, ppn.SV_14, ppn.SV_30, ppn.SV_180, ppn.SV_365]./delayed_reward; % add the delayed reward at the start to calculate the are of the first trapeziod
            diff_delay_norm = diff(delay_norm); % (x2-x1)
            for i=1:length(SV_norm)-1
                sum_SV_norm(i) = SV_norm(i)+SV_norm(i+1); % (y1+y2)
            end
            sum_SV_norm = sum_SV_norm./2; % (y1+y2)/2
            areas = diff_delay_norm.*sum_SV_norm; % (x2-x1)*(y1+y2)/2
            for i=1:length(areas)
                evalc(sprintf('ppn.AUC_area%d = areas(%d);', i,i));
            end
            ppn.AUC = sum(areas);

            clear SV_norm diff_delay_norm sum_SV_norm areas
            curtab = struct2table(ppn);
            evalc(sprintf('Data.TD.session%d(%d,:) = curtab;',s,p));
            clear cur_data delays now_later values 
        end
        if all(stims(p,:) == [1 0]) % session 1 = real, session 2 = sham
            Data.TD.real(p,:) = Data.TD.session1(p,:);
            Data.TD.sham(p,:) = Data.TD.session2(p,:);
        elseif all(stims(p,:) == [0 1]) % session 2 = real, session 1 = sham
            Data.TD.real(p,:) = Data.TD.session2(p,:);
            Data.TD.sham(p,:) = Data.TD.session1(p,:);
        end
        clear ppn 
        evalc(sprintf('Data.TD.session1.bonus(%d) = bonus;',p));
        evalc(sprintf('Data.TD.session2.bonus(%d) = bonus;',p));
        %clear bonus
    end
    clear p s delay* 
    %save the data
    save(outputfile,'Data')
    outputfile2 = strcat(outputfile,'_',date);
    save(outputfile2,'Data')
    fprintf('DATA SAVED')

    %% GoNoGO task

    % loop over participants
    ppns = fieldnames(Data.GoNoGo);
    ppns=ppns(contains(ppns,'pp'));
    stimulus_time=250;
    prac_trls = 20;
    for p=1:length(ppns)
        % loop over the sessions
        sess = eval(sprintf('fieldnames(Data.GoNoGo.%s)',ppns{p}));
        for s=1:length(sess)
            % store the currently relevant data in a temp variable
            cur_data = eval(sprintf('Data.GoNoGo.%s.%s',ppns{p},sess{s}));
            % get the variables of interest for this dataset
            ppn.Subject         = cur_data.Subject(1);
            ppn.Session         = cur_data.Session(1);
            ppn.Age             = cur_data.Age(1);
            ppn.Sex             = cur_data.Sex(1);
            ppn.trial           = cur_data.trial(:);
            ppn.Cijfer          = cur_data.Cijfer(:);
            ppn.CijferExpACC    = cur_data.CijferExpACC(:); % accurate response during stimulus presentation
            ppn.CijferExpRT     = cur_data.CijferExpRT(:); % RT of accurate response during stimulus presentation
            ppn.isi2ACC         = cur_data.isi2ACC(:); % accurate response but after stimulus presentation
            ppn.isi2RT          = cur_data.isi2RT(:); % RT accurate response but after stimulus presentation
            ppn.CijferExpRESP   = cellstr(cur_data.CijferExpRESP); % accurate response during stimulus presentation
            ppn.isi2RESP        = cellstr(cur_data.isi2RESP); % accurate response during stimulus presentation

            % RESP
            for i=prac_trls+1:length(ppn.isi2RT)-1 % exclude practice trials and the last row
                if strcmp(ppn.CijferExpRESP{i}, '{SPACE}') % check if a response is made during stim presentation
                    ppn.RESP(i,:) = 1;
                elseif strcmp(ppn.isi2RESP{i}, '{SPACE}') % check if a response is made during ISI
                    ppn.RESP(i,:) = 1;
                else % no response
                    ppn.RESP(i,:) = 0;
                end
            end
            
            % ACC
            for i=prac_trls+1:length(ppn.isi2RT)-1 % exclude practice trials and the last row
                if ppn.Cijfer(i)==3 && ppn.RESP(i)==0 
                    ppn.ACC(i,:) = 11; % Correct no Go trial
                elseif ppn.Cijfer(i)==3 && ppn.RESP(i)==1 
                    ppn.ACC(i,:) = 12; % Incorrect no Go trial
                elseif ppn.Cijfer(i)~=3 && ppn.RESP(i)==1 
                    ppn.ACC(i,:) = 22; % Correct Go trial
                elseif ppn.Cijfer(i)~=3 && ppn.RESP(i)==0 
                    ppn.ACC(i,:) = 21; % Incorrect Go trial
                end
            end
            ppn.NoGo_corr = length(ppn.ACC(ppn.ACC==11)); % trial counts
            ppn.NoGo_incorr = length(ppn.ACC(ppn.ACC==12)); % trial counts
            ppn.Go_corr = length(ppn.ACC(ppn.ACC==22)); % trial counts
            ppn.Go_incorr = length(ppn.ACC(ppn.ACC==21)); % trial counts
            
            ppn.NoGo_corr_prop = ppn.NoGo_corr/25;
            ppn.NoGo_incorr_prop = ppn.NoGo_incorr/25;
            ppn.Go_corr_prop = ppn.Go_corr/25;
            ppn.Go_incorr_prop = ppn.Go_incorr/25;
            
            % RT
            % if a response is made during the ISI, make the stim RT the stim duration.
            ppn.RT = ppn.CijferExpRT;
            for i=prac_trls+1:length(ppn.isi2RT)-1 % exclude practice trials and the last row
                if ppn.isi2RT(i)>0 && ppn.CijferExpRT(i)==0 % check if a response is made after stim presentation only
                    ppn.RT(i) = stimulus_time; % add stim dur
                    ppn.RT(i) = ppn.RT(i) + ppn.isi2RT(i); % add ISI RT
                elseif ppn.isi2RT(i)==0 && ppn.CijferExpRT(i)==0 && strcmp(ppn.isi2RESP(i),'{SPACE}') % handle responses 'between' stimulus onset and ISI
                    ppn.RT(i) = stimulus_time; % add stim dur
                end
            end
            ppn.RT(ppn.RT==0)=NaN;
            ppn.meanRT_all = nanmean(ppn.RT);
            ppn.meanRT_Go_corr = nanmean(ppn.RT(ppn.ACC==22));
            ppn.meanRT_NoGo_incorr = nanmean(ppn.RT(ppn.ACC==12));
            
            ppn=rmfield(ppn,{'trial','Cijfer','CijferExpACC','CijferExpRT', 'CijferExpRESP', 'isi2ACC', 'isi2RT','isi2RESP', 'RESP', 'ACC', 'ACC2', 'RT'});
            curtab = struct2table(ppn);
            evalc(sprintf('Data.GoNoGo.session%d(%d,:) = curtab;',s,p));
            clear cur_data 
        end
        if all(stims(p,:) == [1 0]) % session 1 = real, session 2 = sham
            Data.GoNoGo.real(p,:) = Data.GoNoGo.session1(p,:);
            Data.GoNoGo.sham(p,:) = Data.GoNoGo.session2(p,:);
        elseif all(stims(p,:) == [0 1]) % session 2 = real, session 1 = sham
            Data.GoNoGo.real(p,:) = Data.GoNoGo.session2(p,:);
            Data.GoNoGo.sham(p,:) = Data.GoNoGo.session1(p,:);
        end
    end
    
    % save the data
    save(outputfile,'Data')
    outputfile2 = strcat(outputfile,'_',date);
    save(outputfile2,'Data')
    fprintf('DATA SAVED')

    close all

end

% save the data
save(outputfile,'Data')
outputfile2 = strcat(outputfile,'_',date);
save(outputfile2,'Data')
fprintf('DATA SAVED')