% DEFAULT CALLBACK FUNCTION --
% CAN BE USED FOR EVERYTHING THAT IS NOT YET DONE
function not_yet_done(varargin)
    uiwait(warndlg('Sorry, this function has not yet been implemented.',...
        'Oops, not yet done...', 'modal')); uiresume;
end