function [badfiles, readout,counter_good,counter_bad,othergoodtest,othergoodtype,otherbadtest,counter_more_good,previous_good] = analyse_line_profile(type,test_line, badfiles, readout,files,s,counter_good,counter_bad,othergoodtest,othergoodtype,otherbadtest,counter_more_good,previous_good)
test = mean(test_line);

manual_override=0;
%Initial boundaries for defining good files.
if test>0.014 && test<0.025
    think = 'Thorium';
    good = strcmp(type,think);
elseif test>0.035 && test <0.04
    think = 'White L';
    good = strcmp(type,think);
elseif test >0.065 && test <0.08
    think = 'Stellar';
    good = strcmp(type,think);
    counter_good=counter_good+1;
else
    %Here are tests based on previous good identifications
    if previous_good>1
        flag=0;
        for n=1:previous_good
            eval(['lowerlimit=0.9*othergoodtest.v' num2str(n) ';'])
            eval(['upperlimit=1.1*othergoodtest.v' num2str(n) ';'])
            if test > lowerlimit && test < upperlimit && flag==0;
                eval(['think = othergoodtype.v' num2str(n) ';'])
                good = strcmp(type,think);
                flag=1;
                counter_more_good=counter_more_good+1;
            end
        end
    end

    %Here are tests based on previous bad identifications
    if counter_bad>1
        for m=1:counter_bad
            eval(['lowerlimit=0.9*otherbadtest.v' num2str(m) ';'])
            eval(['upperlimit=1.1*otherbadtest.v' num2str(m) ';'])
            if test > lowerlimit && test < upperlimit
                good = 0;
            end
        end
    end
    
    if ~exist('good','var')
    %identify program
    manual_override=0;
    frame_type=type;
    save('test_line.mat','test_line','frame_type')
    uiwait(untitled2)
    load('good.mat');
    if good==1
        previous_good=previous_good+1;
        manual_override=1;
        counter_good=counter_good+1;
        eval(['othergoodtest.v' num2str(previous_good) '=test;'])
        eval(['othergoodtype.v' num2str(previous_good) '=type;'])
    elseif good==0
        counter_bad=counter_bad+1;
        eval(['otherbadtest.v' num2str(counter_bad) '=test;'])
    end
end


%check reverse image by looking at max
part_line=find(test_line==1);
if part_line(1)<2032 && manual_override==0
    good = 0;
end

%compares whether it is labelled right.
if good == 0
    badfiles =  cat(1,badfiles, [{strcat(cd,'/', files.list{s}, '.fit')}]);
    readout = cat(1,readout, [{strcat(files.list{s}, ' a  ', type, ' is bad')}]);
    save(strcat('bad_', type, files.list{s}), 'test_line')
end

if good == 1
    save(strcat('good_', type, files.list{s}), 'test_line')
    readout = cat(1,readout, [{strcat(files.list{s}, ' a  ', type, ' is good')}]);
end

end


