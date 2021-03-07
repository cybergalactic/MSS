function vessel = read_veres_ABC(filename, disp_flag)
% vessel = read_veres_ABC(filename, disp_flag)
%
% Read data from ShipX (Veres) output file *.re7
%
% Inputs:
%    filename:                 Text string with name of VERES output file
%    disp_flag (optionally):   1 for displaying VERES run info to screen
%			                   0 turns off display
% Outputs:
% vessel contains data in the following form:
%
%    vessel.main.
%       m:     mass
%       rho    density of water
%       nabla: displacement volume
%       k44:   radius of inertia 
%       k55:   radius of inertia 
%       k66:   radius of inertia 
%       m:     mass
%       CG:    centre of gravity  [LCG 0 VCG] w.r.t. Lpp/2 and Keel Line
%       Lpp:   length between the perpendiculars
%       Lwl:   length of water line
%       T:     draft (water line)    
%       B:     breadth of ship
%          
% Author:    T. I. Fossen
% Date:      2004-06-24 
% Revisions: 2008-02-15 R1.0
%            2013-07-09 Fixed bug, wrong coordinate origin fixed
%            2021-03-07 Bug fixes

if ~exist('disp_flag')
	disp_flag = 0;
end

Tmtrx = diag([-1 1 -1 -1 1 -1]);  % veres2fossen axes

%% ------------------------------------------------------------------------
% Read data from file
%--------------------------------------------------------------------------

% Open file
fid = fopen(filename);

if fid <= 0
	disp('Invalid filename')
end

% Read the first 6 lines, which contain written run info
for m = 1:6	
	header{m} = fgetl(fid);
end

% Read the next 5 lines, which contain numeric run info
for m = 7:11
	header{m} = fgets(fid);
end

% check if this is version 1.0 or 2.0
temp = str2num(header{11});
[dimx,dimy] = size(temp);
if dimy == 6
    version = 1;
    disp('Veres data file *.re7 is version 1.0, processing data...')    
else
    version = 2;
    disp('Veres data file *.re7 is version 2.0, processing data...')    
end

% Line 7: rho and g
temp = str2num(header{7});
rhow = temp(1);
g = temp(2);

% Line 8: Lpp, B, T
temp = str2num(header{8});
LPP = temp(1);
B = temp(2);
T_WL = temp(3);

vessel.main.g   = g;
vessel.main.Lpp = LPP;
vessel.main.T   = T_WL;
vessel.main.B   = B;
vessel.main.rho = rhow;

% Line 9: LCG, VCG
temp = str2num(header{9});
LCG =  temp(1); 
VCG =  temp(2);

vessel.main.CG = [-LCG 0 VCG];  % x postive forwards

if version == 2
    % Line 10: new in version 2
    temp = str2num(header{10});

    % Line 11: nvel, nhead, nfreq, ndof
    temp = str2num(header{11});
    nvel = temp(1);
    nhead = temp(2);
    nfreq = temp(3);
    ndof = temp(4);

    % Read the next 6 lines, which contain rigid body inertia matrix MRB
    for m = 1:6
        mdata{m} = fgets(fid);
    end
else
    % Line 11: nvel, nhead, nfreq, ndof
    temp = str2num(header{10});
    nvel = temp(1);
    nhead = temp(2);
    nfreq = temp(3);
    ndof = temp(4);

    % Read the next 5 lines, which contain rigid body inertia matrix MRB
    mdata{1} = header{11};  % first line of MRB
    for m = 1:5
        mdata{m+1} = fgets(fid);
    end
end

% Extract MRB (computed in waterline just above CG)
MRB = [ str2num(mdata{1})
        str2num(mdata{2})
        str2num(mdata{3})
        str2num(mdata{4})
        str2num(mdata{5})
        str2num(mdata{6})];
    
MRB = Tmtrx*MRB*Tmtrx;         % MRB transformed from veres to fossen axes
H = Hmtrx([-LCG,0,0]);         % Transform MRB to CO midships
vessel.MRB = H'*MRB*H;

vessel.main.m = MRB(1,1);
vessel.main.k44 = sqrt(MRB(4,4)/MRB(1,1));
vessel.main.k55 = sqrt(MRB(5,5)/MRB(1,1));
vessel.main.k66 = sqrt(MRB(6,6)/MRB(1,1));
vessel.main.nabla = MRB(1,1)/rhow;

% Loop over velocities
for velno = 1:nvel

    temp = str2num(fgets(fid));
    vels(velno) = temp(1);

    % Loop over headings
    for headno = 1:nhead

        headings(headno) = str2num(fgets(fid));

        % Loop over frequencies
        for freqno = 1:nfreq

            freqs(freqno) = str2num(fgets(fid));

            Atemp = [str2num(fgets(fid))
                str2num(fgets(fid))
                str2num(fgets(fid))
                str2num(fgets(fid))
                str2num(fgets(fid))
                str2num(fgets(fid))];

            if version == 2
                AtempADD = [str2num(fgets(fid))
                    str2num(fgets(fid))
                    str2num(fgets(fid))
                    str2num(fgets(fid))
                    str2num(fgets(fid))
                    str2num(fgets(fid))];
            end
            
            Btemp = [str2num(fgets(fid))
                str2num(fgets(fid))
                str2num(fgets(fid))
                str2num(fgets(fid))
                str2num(fgets(fid))
                str2num(fgets(fid))];

            if version == 2
                BtempADD = [str2num(fgets(fid))
                    str2num(fgets(fid))
                    str2num(fgets(fid))
                    str2num(fgets(fid))
                    str2num(fgets(fid))
                    str2num(fgets(fid))];
            end
            
            Ctemp = [str2num(fgets(fid))
                str2num(fgets(fid))
                str2num(fgets(fid))
                str2num(fgets(fid))
                str2num(fgets(fid))
                str2num(fgets(fid))];
            
            if version == 2
                CtempADD = [str2num(fgets(fid))
                    str2num(fgets(fid))
                    str2num(fgets(fid))
                    str2num(fgets(fid))
                    str2num(fgets(fid))
                    str2num(fgets(fid))];
            end
                        
            Rolltemp=str2num(fgets(fid));

            Amtrx(:,:,freqno,headno, velno) = Atemp;
            Bmtrx(:,:,freqno,headno, velno) = Btemp;
            Cmtrx(:,:,freqno,headno, velno) = Ctemp;
            Rollvect(:,:,freqno,velno) = Rolltemp;
   
		end % Loop over frequencies
		
	end % Loop over headings
	
end % Loop over velocities

headno = findstr(headings,90);  % 90 deg for data such that w_o=w_e
if isempty(headno), headno = 1; end

% remove heading numbers since A,B,C are independent of heading
Amtrx2 = reshape(Amtrx(:,:,:,headno,:),6,6,nfreq,velno); 
Bmtrx2 = reshape(Bmtrx(:,:,:,headno,:),6,6,nfreq,velno);
Cmtrx2 = reshape(Cmtrx(:,:,:,headno,:),6,6,nfreq,velno);

% transform to Fossen axes
for k = 1:nvel
    for i = 1:nfreq
        vessel.A(:,:,i,k) = H'*Tmtrx*Amtrx2(:,:,i,k)*Tmtrx*H;  
        vessel.B(:,:,i,k) = H'*Tmtrx*Bmtrx2(:,:,i,k)*Tmtrx*H;  
        vessel.C(:,:,i,k) = H'*Tmtrx*Cmtrx2(:,:,i,k)*Tmtrx*H; 
        
        % include viscous roll damping
        Bv1  = Rollvect(1,1,i,k); % linear damping
        Bv2L = Rollvect(1,3,i,k); % nonlinear damping (linearized)
        
        if Bv1 < 0     % remove negative damping for V-shaped hulls
            Bv1 = 0;   % limitation of IKEDA theory
        end
        if Bv2L < 0    % remove negative damping
            Bv2L = 0;
        end
            
        vessel.roll.Bv44(i,k)  = Bv1 + Bv2L;
        
    end
end

vessel.roll.veres   = Rollvect;
vessel.freqs        = freqs;
vessel.headings     = [0:10:350]*pi/180;
vessel.velocities   = vels;

%% ------------------------------------------------------------------------
% ADDED MASS and DAMPING TRANSFORMATIONS
% adds A11 and B11 terms for Veres
% computes viscous friction Bv
%--------------------------------------------------------------------------
[Anew,Bnew,Bv]   = ABCtransform(vessel,'veres',disp_flag);
vessel.Bv        = Bv; 
vessel.A         = Anew;
vessel.B         = Bnew;

%% ------------------------------------------------------------------------
% Display run info on screen
%--------------------------------------------------------------------------
if disp_flag>0
    
    disp('*******************************************************************')
    disp('SHIPX (VERES) DATA')
    disp('*******************************************************************')
    disp(header{2})
    disp(header{3})
    disp(header{4})
    disp(header{5})
    disp(header{6})
    disp('-------------------------------------------------------------------')
    disp(' Run info:')
    disp(' Main particulars:')
    fprintf(' Lpp    : %0.2f m\n',LPP)
    fprintf(' Breadth: %0.2f m\n',B)
    fprintf(' Draught: %0.2f m\n',T_WL)
    disp(' ')
    fprintf(' Number of velocities : %d\n',nvel)
    fprintf(' Number of headings   : %d\n',nhead)
    fprintf(' Number of frequencies: %d\n',nfreq)
    disp('-------------------------------------------------------------------')
    disp([' Headings (deg)     : ' sprintf('%1.0f ',headings)])
    disp([' Frequencies (rad/s): ' sprintf('%1.1f ',freqs)])
    disp([' Velocities (m/s)   : ' sprintf('%1.1f ',vels)])
    disp('-------------------------------------------------------------------')

end