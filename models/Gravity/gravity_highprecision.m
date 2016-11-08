function [acceleration, Jacobian] = gravity_highprecision(p, maxdeg, model) %#eml



%  GRAVITYSPHERICALHARMONIC Implement a spherical harmonic representation
%   of planetary gravity.
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P ) implements the mathematical
%   representation of spherical harmonic planetary gravity based on
%   planetary gravitational potential. Using P, a M-by-3 array of
%   Planet-Centered Planet-Fixed coordinates, GX, GY and GZ, arrays of M
%   gravity values in the x-axis, y-axis and z-axis of the Planet-Centered
%   Planet-Fixed coordinates are calculated for planet using 120th degree
%   and order spherical coefficients for EGM2008 by default.
%
%   Alternate formats for calling zonal harmonic gravity are:
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P, DEGREE )
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P, MODEL )
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P, MODEL, DEGREE )
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P, MODEL, DEGREE, ACTION )
%   [GX GY GZ] = GRAVITYSPHERICALHARMONIC( P, 'Custom', DEGREE, {DATAFILE DFREADER}, ACTION )
%
%   Inputs for spherical harmonic gravity are:
%   P        :a M-by-3 array of Planet-Centered Planet-Fixed coordinates in
%            meters where the z-axis is positive towards the North Pole. For
%            Earth this would be ECEF coordinates.
%   MODEL    :a string specifying the planetary model:
%            'EGM2008' (Earth), 'EGM96' (Earth), 'LP100K' (Moon), 'LP165P'
%            (Moon), 'GMM2B' (Mars), or 'Custom'.  The default is 'EGM2008'.
%   DEGREE   :a scalar value specifying the degree and order of the
%            harmonic gravity model. For 'EGM2008', the maximum degree and
%            order is 2159 and the default degree and order is 120.   For
%            'EGM96', the maximum degree and order is 360 and the default
%            degree and order is 70.  For 'LP100K', the maximum degree and
%            order is 100 and the default degree and order is 60.  For
%            'LP165P', the maximum degree and order is 165 and the default
%            degree and order is 60.  For 'GMM2B', the maximum degree and
%            order is 80 and the default degree and order is 60.  For
%            'Custom', the default degree and order is the maximum degree.
%   DATAFILE :a file containing the planetary gravitational parameter,
%            planet equatorial radius, maximum degree, and normalized
%            spherical harmonic coefficient matrices.
%   DFREADER :a function handle to an MATLAB(R) function which reads
%            DATAFILE.  The MATLAB function must output planetary
%            gravitational parameter in meters cubed per second squared,
%            planet equatorial radius in meters, maximum degree, and the
%            normalized spherical harmonic coefficient matrices, C and S.
%   ACTION   :a string to determine action for out of range input. Specify
%            if out of range input invokes a 'Warning', 'Error', or no
%            action ('None'). The default is 'Warning'.
%
%   Output calculated for the spherical harmonic gravity includes:
%   GX     :an array of M gravity values in the x-axis of the
%          Planet-Centered Planet-Fixed coordinates in meters per second
%          squared.
%   GY     :an array of M gravity values in the y-axis of the
%          Planet-Centered Planet-Fixed coordinates in meters per second
%          squared.
%   GZ     :an array of M gravity values in the z-axis of the
%          Planet-Centered Planet-Fixed coordinates in meters per second
%          squared.
%
%   Limitations:
%
%   This function has the limitations of excluding the centrifugal effects
%   of planetary rotation, and the effects of a precessing reference frame.
%
%   Spherical harmonic gravity model is valid for radial positions greater
%   than the planet equatorial radius.  Using it near or at the planetary
%   surface can probably be done with negligible error.  The spherical
%   harmonic gravity model is not valid for radial positions less than
%   planetary surface.
%
%   Examples:
%
%   Calculate the gravity in the x-axis at the equator on the surface of
%   Earth, using the 120 degree model of EGM2008 with warning actions:
%       gx = gravitysphericalharmonic( [-6378.1363e3 0 0] )
%
%   Calculate the gravity at 25000 meters over the south pole of Earth using
%   the 70 degree model of EGM96 with error actions:
%       [gx, gy, gz] = gravitysphericalharmonic( [0 0 -6381.751e3], 'EGM96', 'Error' )
%
%   Calculate the gravity at 15000 meters over the equator and 11000 meters
%   over the north pole using a 30th order GMM2B Mars model with warning
%   actions:
%       p  = [2412.648e3 -2412.648e3 0; 0 0 3376.2e3]
%       [gx, gy, gz] = gravitysphericalharmonic( p, 'GMM2B', 30, 'Warning' )
%
%   Calculate the gravity at 15000 meters over the equator and 11000 meters
%   over the north pole using a 60th degree custom planetary model with no
%   actions:
%       p       = [2412.648e3 -2412.648e3 0; 0 0 3376e3]
%       [gx, gy, gz] = gravitysphericalharmonic( p, 'custom', 60, ...
%                       {'GMM2BC80_SHA.txt' @astReadSHAFile}, 'None' )
%
%   See also GRAVITYWGS84, GRAVITYCENTRIFUGAL, GRAVITYZONAL, GEOIDEGM96

%   References:
%   [1] Vallado, D. A., "Fundamentals of Astrodynamics and Applications",
%       McGraw-Hill, New York, 1997.
%   [2] NIMA TR8350.2: "Department of Defense World Geodetic System 1984,
%       Its Definition and Relationship with Local Geodetic Systems."
%   [3] Konopliv, A. S., S. W. Asmar, E. Carranza, W. L. Sjogen, D. N.
%       Yuan., "Recent Gravity Models as a Result of the Lunar Prospector
%       Mission", Icarus, Vol. 150, no. 1, pp 1?18, 2001.
%   [4] Lemoine, F. G., D. E. Smith, D.D. Rowlands, M.T. Zuber, G. A.
%       Neumann, and D. S. Chinn, "An improved solution of the gravity
%       field of Mars (GMM-2B) from Mars Global Surveyor", J. Geophys. Res.,
%       Vol. 106, No. E10, pp 23359-23376, October 25, 2001.
%   [5] Kenyon S., J. Factor, N. Pavlis, and S. Holmes, "Towards the Next
%       Earth Gravitational Model", Society of Exploration Geophysicists
%       77th Annual Meeting, San Antonio, Texas, September 23-28, 2007.
%   [6] Pavlis, N.K., S.A. Holmes, S.C. Kenyon, and J.K. Factor, "An Earth
%       Gravitational Model to Degree 2160: EGM2008", presented at the 2008
%       General Assembly of the European Geosciences Union, Vienna,
%       Austria, April 13-18, 2008.
%   [7] Grueber, T., and A. Kohl, "Validation of the EGM2008 Gravity Field
%       with GPS-Leveling and Oceanographic Analyses", presented at the IAG
%       International Symposium on Gravity, Geoid & Earth Observation 2008,
%       Chania, Greece, June 23-27, 2008.

    %% Initialize

    persistent GM Re degree C S %#ok
    persistent prevmodel
    if isempty(prevmodel)
        prevmodel = 0; end 

    switch model 

        % [GM, Re, degree, C, S]
        case 1%'egm2008'
            % Earth            
            if prevmodel~=1
                load('aeroegm2008.mat')                 
                prevmodel = 1;
            end
            
        case {0 2}%'egm96'
            % Earth            
            if prevmodel~=2
                load('aeroegm96.mat') 
                prevmodel = 2;
            end
          
            %{
        case 3%'lp100k'
            % Moon            
            load('aerolp100k.mat') 
            
        case 4%'lp165p'
            % Moon            
            load('aerolp165p.mat') 
            
        case 5%'gmm2b'
            % Mars            
            load('aerogmm2b.mat') 
            %}
            
    end

    % convert to spherical coordinates
    r      = norm(p);
    phi    = asin (p(3)/r);
    lambda = atan2(p(2),p(1));

    % Check if we're flying underground
    if r <= Re
        % TODO
    end
    
    % Often-used values
    rm1 = 1/r;      pio2  = pi/2;             
    rm2 = rm1*rm1;  sqrt3 = sqrt(3);

    
    %% Fully normalized associated Legendre polynomials 
    
    % Computed via recursion relations
    
    P = zeros(maxdeg+3, maxdeg+3);
    
    cphi = cos(pio2 - phi);    cphi(abs(cphi) <= eps) = 0;
    sphi = sin(pio2 - phi);    sphi(abs(sphi) <= eps) = 0;
    
    P(1,1) = 1;            
    P(2,1) = sqrt3 * cphi; 
    P(2,2) = sqrt3 * sphi; 

    for n = 2 : maxdeg+2
        
        k  = n + 1;
        n2 = 2*n;
        
        sn2   = sqrt(n2+0);  sn2p1 = sqrt(n2+1); 
        sn2m3 = sqrt(n2-3);  sn2m1 = sqrt(n2-1);
                
        P(k,1) = sqrt(n2+1)/n * ...
                 (sn2m1*cphi*P(k-1,1) - (n-1)/sn2m3 * P(k-2,1));             
        P(k,k) = sn2p1/sn2 * sphi*P(k-1,k-1);  
        
        for m = 1 : n-1
            q = m + 1;            
            P(k,q) = sn2p1 / (sqrt(n+m)*sqrt(n-m)) * ...
                     (sn2m1*cphi*P(k-1,q) - sqrt(n+m-1)*sqrt(n-m-1)/sn2m3*P(k-2,q));
        end
        
    end
    
    
    %% Compute partial derivatives of potential 
    
    % pre-compute trig factors (recursion relations)
    slambda  = sin(lambda);  
    clambda  = cos(lambda);  
    
    smlambda = [slambda  2*clambda*slambda     zeros(1, maxdeg-1)];    
    cmlambda = [clambda  2*clambda*clambda-1   zeros(1, maxdeg-1)];
    
    for m = 3:maxdeg+1
        smlambda(m) = 2*clambda*smlambda(m-1) - smlambda(m-2);
        cmlambda(m) = 2*clambda*cmlambda(m-1) - cmlambda(m-2);
    end
    
    % Other often-used values    
    roR  = Re*rm1;             roR_n = roR;
    xy2  = p(1)^2 + p(2)^2;    nxy   = sqrt(xy2);
    oxy2 = 1/xy2;              onxy  = 1/nxy;

    % Summations in radial coordinates    
    dUdr_n      = 1;
    dUdphi_n    = 0;
    dUdlambda_n = 0;
    for n = 2:maxdeg
        
        k     = n+1;
        roR_n = roR_n * roR;
        
        % (summation over m)
        fac1        = ( C(k,1:k).*cmlambda(1:k) + S(k,1:k).*smlambda(1:k) ).';
        fac2        = ( S(k,1:k).*cmlambda(1:k) - C(k,1:k).*smlambda(1:k) ).';
        dUdr_m      = P(k,1:k) * fac1;
        dUdphi_m    = (P(k,2:k+1) - p(3)*onxy * (0:n).*P(k,1:k)) * fac1;
        dUdlambda_m = (0:n).*P(k,1:k) * fac2;
        
        % (summation over n)
        dUdr_n      = dUdr_n      + k*roR_n*dUdr_m;     
        dUdphi_n    = dUdphi_n    +   roR_n*dUdphi_m;   
        dUdlambda_n = dUdlambda_n +   roR_n*dUdlambda_m;
        
    end

    % Acceleration in spherical coordinates
    dUdr      = -GM*rm2 * dUdr_n;
    dUdphi    =  GM*rm1 * dUdphi_n;
    dUdlambda =  GM*rm1 * dUdlambda_n;
    
    % Special case for positions near the poles
    if abs(atan2(p(3),nxy)) == pio2
        
        acceleration = [
            0
            0
            p(3)*dUdr*rm1;
            ];
    
    % All other locations:
    else   
        % Convert back to ECEF
        f = dUdr*rm1;
        g = f - p(3)*dUdphi*rm2*onxy;
        h = dUdlambda*oxy2;
        
        acceleration = [
            g*p(1) - h*p(2)
            g*p(2) + h*p(1)
            f*p(3) + nxy*dUdphi*rm2
            ];
    end

    
    %% Compute partial derivatives (Jacobian) 

end



