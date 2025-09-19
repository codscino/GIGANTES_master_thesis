function [r1dot,r2dot]=LambertArc_IntraNode_Transfer(r1, r2, r1guess, TOF, mu)




%% Step 2: Differential Corrector
% The algorithm will continuate a time vetor until matching the expected
% arrival time.
ToF_target=TOF; % Required arrival time; target transfer time.

% Seed orbit
sma_seed=1e9;
r1dot_seed=r1guess;

%--------------------------------------------------------------------------
% Defining the number of differential correction steps
nIterations=25; % number of Iterations until reaching the desired solution
% Start loop

    DT=ToF_target;
    %--------------------------------------------------------------------------
    % STEP 4: Shooting method
    ETolerance=1e-3; % 1 m miss match are allowed
    numMaxIter=25;
    numIter=0;
    r2_Iteration=0;
    while (norm(r2_Iteration-r2)>ETolerance)&&(numIter<numMaxIter)
        
        % [r2Star]=FGKepler_dt (r1, r1dot, DT)
        % [Phi_12]=STM_Lambert (r1, r1dot_seed, DT)
        
        [r2Star, Phi_12]=FGKepler_dt_STM_Lambert (r1, r1dot_seed, DT, mu);
        r2_Iteration=r2Star;
        Dr1dot=Phi_12\(r2-r2Star)';
        
        r1dot_seed=r1dot_seed+Dr1dot';
        
        numIter=numIter+1;
        
        % checks and warnings;
        sma_seed=mu/2/(mu/norm(r1)-norm(r1dot_seed)^2/2);
        if sma_seed<0
            warning('sma_seed error:The shooting method did not converge ')
            r1dot=[NaN NaN NaN];
            r2dot=[NaN NaN NaN];
            return
        end
    end
    if numIter>=numMaxIter
            warning('numMaxIter issue:The shooting method did not converge')
            r1dot=[NaN NaN NaN];
            r2dot=[NaN NaN NaN];
            return
    end

%--------------------------------------------------------------------------
%% Step 5:End of Algorithm
% how much is Df for the corrsponding DT
sma_seed=mu/2/(mu/norm(r1)-norm(r1dot_seed)^2/2);
n=sqrt(mu/sma_seed^3);
DM=n*DT;
sigma0=dot(r1,r1dot_seed)/sqrt(mu);
fun=@(DE)DM-DE-sigma0/sqrt(sma_seed)*(1-cos(DE))+(1-norm(r1)/sma_seed)*sin(DE);
DE_f = fzero(fun,DM);

% defining Transition matrix
% F
F=1-sma_seed/norm(r1)*(1-cos(DE_f));
% F_old=1-r2_n/p_min*(1-cos(Df));
% G
G=DT+sqrt(sma_seed^3/mu)*(sin(DE_f)-DE_f);
% G_old=r1_n*r2_n/sqrt(MuSun*p_seed)*sin(Df);
% Gdot
rf=sma_seed+(norm(r1)-sma_seed)*cos(DE_f)+sqrt(sma_seed)*sigma0*sin(DE_f);
Gdot=1-sma_seed/rf*(1-cos(DE_f));
%--------------------------------------------------------------------------
% STEP 4: Shooting method

% propagate from r1 to r2 in the step time corrected
r1dot=r1dot_seed;
r2=F*r1+G*r1dot; %watch out since you are modifying things
r2dot=1/G*(-r1+Gdot*r2);


end
%% Auxiliary function minETransfer
function [sma_minE, ToF_min, r1dot]=minETransfer(r1, r2, tm, mu)

%--------------------------------------------------------------------------
% Semimajor Axis of minimum energy orbit transfer
r1_n=norm(r1);
r2_n=norm(r2);

c=norm(r2-r1); %chord
sma_minE=1/4*(r1_n+r2_n+c);

% angular distance between the points
cos_DtrA=(r2_n^2+r1_n^2-c^2)/(2*r1_n*r2_n);
sin_DtrA=tm*sqrt(1-cos_DtrA^2);

Delta_trA=acos((r2_n^2+r1_n^2-c^2)/(2*r1_n*r2_n));
if tm<0
    Delta_trA=2*pi-Delta_trA;
end

% Semilatus rectum as Battin pg 246
p_min=r1_n*r2_n/c*(1-cos_DtrA);
% how much is the eccentricity of this orbit?
ecc_min=sqrt(1-p_min/sma_minE);

% how long it takes to go from r1 to r2?

n=sqrt(mu/sma_minE^3);
% ToF_min=DM/n; % Time of Flight in Seconds
% ToF_minAsVallado=1/3*sqrt(2/mu)*((2*sma_minE)^(3/2)-(2*sma_minE-c)^(3/2))
betaE=2*asin(sqrt((2*sma_minE-c)/2/sma_minE));
ToF_min=1/n*(pi-sign(tm)*(betaE-sin(betaE)));
% TimeMissMatch_days=TOF-ToF_min/3600/24

%--------------------------------------------------------------------------
% compute r1dot
% F
F=1-r2_n/p_min*(1-cos_DtrA);
% G
G=r1_n*r2_n/sqrt(mu*p_min)*sin_DtrA;
%r1dot
r1dot=1/G*(r2-F*r1);

end
%% Auxiliary function FGKepler_dt_STM_Lambert
function [r2Star, Phi_12]=FGKepler_dt_STM_Lambert (r1, r1dot_seed, DT, mu)

r1_n=norm(r1);
% vis viva
sma_seed=mu/2/(mu/r1_n-norm(r1dot_seed)^2/2);
% how much is Df for the corrsponding DT
n=sqrt(mu/sma_seed^3);

DM=n*DT;
sigma0=dot(r1,r1dot_seed)/sqrt(mu);
fun=@(DE)DM-DE-sigma0/sqrt(sma_seed)*(1-cos(DE))+(1-norm(r1)/sma_seed)*sin(DE);
DE_f = fzero(fun,DM);

% defining Transition matrix
% F
F=1-sma_seed/r1_n*(1-cos(DE_f));
% F_old=1-r2_n/p_min*(1-cos(Df));
% G
G=DT+sqrt(sma_seed^3/mu)*(sin(DE_f)-DE_f);
%--------------------------------------------------------------------------
% STEP 4: Shooting method

% propagate from r1 to r2 in the step time corrected
r2Star=F*r1+G*r1dot_seed; %watch out since you are modifying things

% Relevant State Transition Matrix (as in page 432 from Schaub)
C=sma_seed*sqrt(sma_seed^3/mu)*(3*sin(DE_f)-(2+cos(DE_f))*DE_f)...
    -DT*sma_seed*(1-cos(DE_f));

dr=r2Star-r1;
%r2dot
% Gdot
Gdot=1-sma_seed/norm(r2Star)*(1-cos(DE_f));

r2dot_iteration=1/G*(-r1+Gdot*r2Star);
dv=r2dot_iteration-r1dot_seed;

Phi_12=r1_n/mu*(1-F)*(dr'*r1dot_seed - dv'*r1)+C/mu*r2dot_iteration'*r1dot_seed+G*eye(3);

end