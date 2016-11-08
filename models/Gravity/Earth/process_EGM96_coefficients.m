
% (This is the old "Init_ENV_EGR")
function process_EGM96_coefficients()

    %% Parameters
    
    maxdegree  = 360 + 1;
    maxdotdeg  = 3;
    
    % (from readme)
    GM         = 3986004.415e8;
    Re         = 6378136.3;
    refJD      = 2446431.5;  % 1 Jan 1986
    
    Cnm_dot    = [              0                0  0 
                                0                0  0 
                  +1.16275534e-11  -0.32000000e-11  0];
    Snm_dot   =  [              0                0  0 
                                0                0  0 
                                0  +1.62000000e-11  0];                      
   
    pth        = fileparts(mfilename('fullpath'));
    inputfile  = fullfile(pth, 'egm96_to360.ascii');
    CFile      = fullfile(pth, 'EGM96_parameters.c');
    HFileName  = 'EGM96_parameters.h';
    HFile      = fullfile(pth, HFileName);
    
    prefix     = 'Gravity_';
    
    GMname     = [prefix 'GM'];
    REname     = [prefix 'Re'];
    
    LegendreP  = [prefix 'P'];
    smlambda   = [prefix 'smlambda'];
    cmlambda   = [prefix 'cmlambda'];
    
    Cnmname    = [prefix 'Cnm'];
    Snmname    = [prefix 'Snm'];
    
    Cnmdotname = [prefix 'Cnm_dot'];
    Snmdotname = [prefix 'Snm_dot'];
    refdate    = [prefix 'reference_epoch'];    
    
    maxdegtxt  = [upper(prefix) 'COEFF_MAXDEGREE'];
    maxddrftxt = [upper(prefix) 'DRIFT_MAXDEGREE'];
    
    
    %% Read and format Cnm/Snm parameters and rates
    
    % Extract and reshape data     
    grav_param = load('-ASCII', inputfile);
        
    datasz = grav_param(end,1)+1;
    n      = grav_param(:,1)+1;
    m      = grav_param(:,2)+1;    
    inds   = sub2ind([datasz datasz], n,m);
    
    Cnm = zeros(datasz);  Cnm (inds) =  grav_param(:,3);
    Snm = zeros(datasz);  Snm (inds) =  grav_param(:,4);
    
    % Re-format to C style array definition
    Cnm = format_coeffs(Cnm);
    Snm = format_coeffs(Snm);
    
    % (also the rates)
    Cnm_dot = format_coeffs(Cnm_dot);
    Snm_dot = format_coeffs(Snm_dot);    
    
    %% Write header file
    
    header_guardtxt = upper(regexprep(HFileName, '\.', '_'));    
    
    content = {
        sprintf('#ifndef %s\n#define %s\n', header_guardtxt, header_guardtxt)        
        ' '
        sprintf('#define %s (%du)', maxdegtxt , maxdegree)
        sprintf('#define %s (%du)', maxddrftxt, maxdotdeg)
        ' '
        sprintf(['double %s[%s+3u][%s+3u],\n',...
                 '       %s[%s+1u],\n',...
                 '       %s[%s+1u];'],...
                 LegendreP, maxdegtxt,maxdegtxt,...
                 cmlambda, maxdegtxt,...
                 smlambda, maxdegtxt);
        ' '     
        sprintf( 'extern const double %s;', GMname);
        sprintf( 'extern const double %s;', REname);
        ' '
        sprintf('extern const double %s[%s][%s];', Cnmname, maxdegtxt, maxdegtxt)
        sprintf('extern const double %s[%s][%s];', Snmname, maxdegtxt, maxdegtxt)
        ' '
        sprintf('extern const double %s;', refdate)
        sprintf('extern const double %s[%s][%s];', Cnmdotname, maxddrftxt, maxddrftxt)
        sprintf('extern const double %s[%s][%s];', Snmdotname, maxddrftxt, maxddrftxt)                
        ' '
        '#endif'
        };
    
    fid = fopen(HFile,'w');
    OC2 = onCleanup(@() any(fopen('all')==fid) && fclose(fid));
    fprintf(fid, '%s\n', content{:});
    fclose(fid);
    
    
    %%  Write C-file
              
    content = {...   
        sprintf('#include "%s"', HFileName)
        ' '
        sprintf( 'const double %s = %f;', GMname, GM);
        sprintf( 'const double %s = %f;', REname, Re);            
        ' '
        sprintf( 'const double %s[%s][%s] = %s', Cnmname, maxdegtxt, maxdegtxt, sprintf('%s\n', Cnm{:}) )            
        sprintf( 'const double %s[%s][%s] = %s', Snmname, maxdegtxt, maxdegtxt, sprintf('%s\n', Snm{:}) )        
        ' '        
        sprintf( 'const double %s = %f;', refdate, refJD)      
        ' '
        sprintf( 'const double %s[%s][%s] = %s', Cnmdotname, maxddrftxt, maxddrftxt, sprintf('%s\n', Cnm_dot{:}) )            
        sprintf( 'const double %s[%s][%s] = %s', Snmdotname, maxddrftxt, maxddrftxt, sprintf('%s\n', Snm_dot{:}) )        
        };
        
    fid = fopen(CFile,'w');
    OC1 = onCleanup(@() any(fopen('all')==fid) && fclose(fid));
    fprintf(fid, '%s\n', content{:});
    fclose(fid);
    
end


function str = format_coeffs(C)

    N = size(C,1);
    str = cell(N,1);
    for ii = 1:N
        str{ii} = ['    {',...
            sprintf('%+.16e, ', C(ii,1:end-1)),...
            sprintf('%+.16e},', C(ii,end))];
    end

    % Strip last comma
    str{end} = str{end}(1:end-1);

    % Finish off C style 2D array 
    str = [
        '{'
        str
        '};'];
    
end

