% counts the total number of lines in a file, very quickly
function numlines = countlines(filename)
    if (isunix) %# Linux, mac
        [dummy, result] = system( ['wc -l ', filename] ); %#ok<ASGLU>
    else % if (ispc) %# Windows, and any other
        result = perl('countlines.pl', filename);        
    end    
    numlines = str2double(result);
end
