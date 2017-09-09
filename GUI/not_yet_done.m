% DEFAULT CALLBACK FUNCTION --
% CAN BE USED FOR EVERYTHING THAT IS NOT YET DONE
function not_yet_done(~,~)
    uiwait(warndlg('Sorry, this function has not yet been implemented.',...
                   'Not yet implemented',...
                   'modal'));
    uiresume;
end