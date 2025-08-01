function find_good_lines_in_synthlist(wavelengths,widths,elem,rang,ratio,minim)
for s=1:numel(wavelengths)
    inds=find(wavelengths > wavelengths(s)-rang & wavelengths < wavelengths(s)+rang);
    if numel(inds) > 0
        testlist=widths(inds);
        if max(testlist) > widths(s)  %is the testline the strongest in the range?
            continue
        else
            test=find(testlist > widths(s)*ratio);
            if numel(test) > 1 | widths(s) < minim
                %are there any other lines greater than the testline*0.1 or
                %is this line too small?
                continue
            else %we have a winner
                fprintf('%.2f A %s line, width=%.2f mA\n',wavelengths(s),elem{s},widths(s))
            end
        end
    end
end

        