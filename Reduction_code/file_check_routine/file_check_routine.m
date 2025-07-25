function [change_flag,counter_good]=file_check_routine(starname)

%%%%%%   PARAMETERS   %%%%%%%
%leave these, used later
badfiles =[];
readout =[];


%%%%%%   MAIN SECTION %%%%%%%
%general bad file routine
has_run=badfile_check_general;
if has_run==0
    info = display_directory_file_information_nooutput();
    save('file_info','info')
else
    load file_info.mat
end

%calibration bad file routine
change_flag=0;
%CYCLE THROUGH FILES FIRST CHECKS FOR CALIBRATIONS AND SKIPS IF NOT
%Get line profile and move it to the right folder
test_badcal=dir('badfiles_calibration.mat');
if numel(test_badcal)>0
    disp('Calibration file check already completed in this directory')
    badfiles_calibration=[];
else
    change_flag=1;
    h = waitbar(0,'Please wait...');
    counter_good=0;
    counter_bad=0;
    counter_more_good=0;
    previous_good=0;
    othergoodtest.v1=[];
    othergoodtype.v1=[];
    otherbadtest.v1=[];
    disp('Testing calibration files')
    badfiles_calibration=[];
    for s=1:numel(info.list)
        waitbar(s/numel(info.list), h, info.list{s})
        %I.D. the file
        [type,star] = identify_the_file(info,s);
        if exist(strcat('god', type, info.list{s},'.mat'),'file')==2
            continue
        end
        if strcmp(type, 'White L') || strcmp(type, 'Thorium')
            [test_line] = get_the_line_profile(star,starname,type,s,info);
            clear 'good'
            [badfiles_calibration, readout,counter_good,counter_bad,othergoodtest,othergoodtype,otherbadtest,counter_more_good,previous_good] = analyse_line_profile(type,test_line, badfiles_calibration, readout,info,s,counter_good,counter_bad,othergoodtest,othergoodtype,otherbadtest,counter_more_good,previous_good);
        end
    end
    char(readout)
    close(h)
    save('badfiles_calibration', 'badfiles_calibration')
end

%CYCLE THROUGH FILES FIRST CHECKS FOR STELLARS AND SKIPS IF NOT
%Get line profile and move it to the right folder
counter_good=0;
counter_more_good=0;
previous_good=0;
counter_bad=0;
badfiles_stellar=[];
%include alias names
wildcard_matching=0;
alt_starnames{1}=starname;
load Mus_known_alias_names
if find(ismember(wildc,starname))>0
    wildcard_matching=1;
    [x,y]=find(ismember(wildc,starname));
    for n=1:length(wildc(x,:))
        alt_starnames{length(alt_starnames)+1}=wildc{x(n),2};
    end
end

for n=1:numel(alt_starnames)
    starname=alt_starnames{n};
    
    badstellar_name=['badfiles_' starname '.mat'];
    test_badstellar=dir(badstellar_name);  
    
    if numel(test_badstellar)>0
        disp(['Stellar file check ' starname ' already completed in this directory'])
        for s=1:numel(info.list)
            %I.D. the file
            [type,star] = identify_the_file(info,s);
            
            if exist(strcat('god', type, info.list{s},'.mat'),'file')==2
                continue
            end
            
            if strcmp(type, 'Stellar') && strcmpi(star,starname)
                counter_good=counter_good+1;
            end
            
        end
    else
        change_flag=1;
        h = waitbar(0,'Please wait...');
        othergoodtest.v1=[];
        othergoodtype.v1=[];
        otherbadtest.v1=[];
        for s=1:numel(info.list)
            waitbar(s/numel(info.list), h, info.list{s})
            %I.D. the file
            [type,star] = identify_the_file(info,s);
            
            if exist(strcat('god', type, info.list{s},'.mat'),'file')==2
                continue
            end
            
            if strcmp(type, 'Stellar') && strcmpi(star,starname)
                [test_line] = get_the_line_profile(star,starname,type,s,info);
                clear 'good'
                [badfiles_stellar, readout,counter_good,counter_bad,othergoodtest,othergoodtype,otherbadtest,counter_more_good,previous_good] = analyse_line_profile(type,test_line, badfiles_stellar, readout,info,s,counter_good,counter_bad,othergoodtest,othergoodtype,otherbadtest,counter_more_good,previous_good);
            end
        end
        char(readout);
        close(h)
        eval(['save(''' badstellar_name ''',''badfiles_stellar'')'])
    end
end

counter_good=counter_good+counter_more_good;

%%%%%%%%%%%%%%%%%%%  MAIN REDUCTION COMPLETE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THE FOLLOWING WILL PLOT EACH OF THE SUB FOLDERS IN LINE ANALYSIS %%%%%%
%%%%%%%%%%%%%%%%%%%%%% SAVE THEM. THIS IS FOR YOU TO INSPECT %%%%%%%%%%%%%%
god_names = dir('good_*.mat');
bad_names = dir('bad_*.mat');

if isempty(bad_names)==0 || isempty(god_names)==0
    figure(1000)
    for s=1:numel(bad_names)
        load(bad_names(s).name);
        sort = bad_names(s).name;
        if sort(5) == 'S'
            subplot(2,3,4)
            plot(test_line)
            hold on
        elseif sort(5) == 'W'
            subplot(2,3,5)
            plot(test_line)
            hold on
        elseif sort(5) == 'T'
            subplot(2,3,6)
            plot(test_line)
            hold on
        end
    end
    
    for s=1:numel(god_names)
        load(god_names(s).name);
        sort = god_names(s).name;
        if sort(6) == 'S'
            subplot(2,3,1)
            plot(test_line)
            hold on
        elseif sort(6) == 'W'
            subplot(2,3,2)
            plot(test_line)
            hold on
        elseif sort(6) == 'T'
            subplot(2,3,3)
            plot(test_line)
            hold on
        end
    end
    
    subplot(2,3,1)
    xlabel(strcat(num2str(numel(dir('good_S*.mat'))),' ',  ' -Stellar Accepted'))
    subplot(2,3,2)
    xlabel(strcat(num2str(numel(dir('good_W*.mat'))),' ',  '-White Lamps Accepted'))
    subplot(2,3,3)
    xlabel(strcat(num2str(numel(dir('good_T*.mat'))),' ',  '-Thoriums Accepted'))
    subplot(2,3,4)
    xlabel(strcat(num2str(numel(dir('bad_S*.mat'))),' ',  '-Rejected'))
    subplot(2,3,5)
    xlabel(strcat(num2str(numel(dir('bad_W*.mat'))),' ',  '-Rejected'))
    subplot(2,3,6)
    xlabel(strcat(num2str(numel(dir('bad_T*.mat'))),' ',  '-Rejected'))
    
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Press enter to move bad files to a subdirectory or CTRL-C to change')
    p=input('break out.');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    try close('1000')
    catch
    end
else
end

try remove_files(badfiles_calibration)
catch
end
try remove_files(badfiles_stellar)
catch
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('The file check program is now complete. The bad files have been')
disp('moved to a separate folder called removed_files and should not be reduced.')
disp('Everything else remaining in the current directory is suitable for further')
disp('use. Reduction will commence shortly.')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end




