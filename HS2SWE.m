function odata = HS2SWE(idata,RhoNew,RhoMax,SnoTemp,Visc,DQMultInc,DQMultMax,HsAcc,c1,c2,c3,c4,c5,g,dt)
%HS2SWE   Convert snow depth recordings to snow water equivalents
%   odata = HS2SWE(idata) where idata is an array with snow depth
%   recordings (cm) with each column representing a time series of snow
%   depth measurements from a station. The output (odata) is an array with
%   simulates snow water equivalents (mm).

arguments
  idata
  RhoNew    =  113.7
  RhoMax    =  571.6
  SnoTemp   = -0.000
  Visc      =  6.051*1e7
  DQMultInc = 0.1
  DQMultMax = 5
  HsAcc     = 2
  c1        = 2.8e-6
  c2        = 0.042
  c3        = 0.046
  c4        = 0.081
  c5        = 0.018
  g         = 9.81
  dt        = 86400
end

PAR.c1    = c1;
PAR.c2    = c2;
PAR.c3    = c3;
PAR.c4    = c4;
PAR.c5    = c5;
PAR.g     = g;
PAR.dt    = dt;

% input data checks/preparation
assert(all(idata(~isnan(idata))>=0),"Negative values in input data")

% looping over time series of data
for serix = 1:size(idata,2)
  disp("Running station: " + serix)
  psdix = find(~isnan(idata(:,serix)) & idata(:,serix) > 0)';
  if isempty(psdix)
  else
    sepix = [0 find(diff(psdix)>1) length(psdix)];
    % looping over continuous strings of non-nan data with HS > 0
    for tsrix = 1:length(sepix)-1
      HS = [0 idata(psdix(sepix(tsrix)+1):psdix(sepix(tsrix+1)),serix)'];

      % initializing data field with a zero height layer @ RhoMax
      MR     = [];
      MR.HS  = zeros(1,length(HS));        % layer depth matrix
      MR.RHO = ones(1,length(HS)).*RhoMax; % layer density matrix
      MR.OVB = zeros(1,length(HS));        % overburden mass matrix
      MR.AGE = 0:length(HS)-1;             % layer age matrix
      MR.DIA = zeros(5,length(HS));        % performance diagnostics

      %loop over days (tn)
      for tn = 2:length(HS)
        %step 1: densification of exisiting layer (limited by RhoMax); ATTN: any changes of the below line need to be replicated in step 3.3
        MR.RHO(:,tn) = min(MR.RHO(:,tn-1)+MR.RHO(:,tn-1).*PAR.dt.*(MR.OVB(:,tn-1).*PAR.g./(Visc.*exp(PAR.c4.*SnoTemp+PAR.c5.*MR.RHO(:,tn-1)))+PAR.c1.*exp(-PAR.c2.*SnoTemp-PAR.c3.*max(0,MR.RHO(:,tn-1)-RhoNew))),RhoMax);
        %step 2: settling of snow according to step 1 assuming constant SWE; %ATTN: any changes of the below line need to be replicated in step 3.3
        MR.HS(:,tn) = MR.HS(:,tn-1)./(MR.RHO(:,tn)./MR.RHO(:,tn-1));
        %step 3: assimilate measured HS (add new snow / melt snow)
        %step 3.0 if HSmeas > HSmod for first time step assume new snow fall and add layer
        if HS(tn) > sum(MR.HS(:,tn)) && tn == 2
          nlix = size(MR.HS,1)+1;
          MR.HS(nlix,:)       = zeros(1,length(HS));
          MR.HS(nlix,tn)      = HS(tn) - sum(MR.HS(:,tn));
          MR.RHO(nlix,:)      = ones(1,length(HS)).*RhoNew;
          MR.OVB(nlix,:)      = zeros(1,length(HS));
          MR.AGE(nlix,:)      = zeros(1,length(HS));
          MR.AGE(nlix,tn:end) = 0:length(HS)-tn;
          if HS(tn) < HS(tn-1)
            MR.DIA(1,tn) = 1; % if observed HS is decreasing while model adds new snow add a note in MR.DIA
          end
          %step 3.1 if HSmeas > HSmod + HSacc assume new snow fall and add layer
        elseif HS(tn) > sum(MR.HS(:,tn)) + HsAcc
          nlix = size(MR.HS,1)+1;
          MR.HS(nlix,:)       = zeros(1,length(HS));
          MR.HS(nlix,tn)      = HS(tn) - sum(MR.HS(:,tn));
          MR.RHO(nlix,:)      = ones(1,length(HS)).*RhoNew;
          MR.OVB(nlix,:)      = zeros(1,length(HS));
          MR.AGE(nlix,:)      = zeros(1,length(HS));
          MR.AGE(nlix,tn:end) = 0:length(HS)-tn;
          if HS(tn) < HS(tn-1)
            MR.DIA(1,tn) = 1; % if observed HS is decreasing while model adds new snow add a note in MR.DIA
          end
          %step 3.2 if HSmeas == HSmod don't do anything
        elseif HS(tn) == sum(MR.HS(:,tn))
          MR.DIA(2,tn) = 0; % note difference between HSmeas - HSmod if positive in MR.DIA
          %step 3.3 if HSmeas > HSmod reapply densification with gradually decreasing densification rate until HSmeas <= HSmod
        elseif HS(tn) > sum(MR.HS(:,tn))
          MR.DIA(2,tn) = HS(tn)-sum(MR.HS(:,tn)); % note difference between HSmeas - HSmod before assimilation
          %step 3.3.1 decreasing deinsification rate
          DQMultCur = 1;
          while mean(MR.RHO(:,tn)) < RhoMax && HS(tn) > sum(MR.HS(:,tn)) && DQMultCur < DQMultMax
            DQMultCur = DQMultCur + DQMultInc;
            MR.RHO(:,tn) = min(MR.RHO(:,tn-1)+MR.RHO(:,tn-1).*PAR.dt.*(MR.OVB(:,tn-1).*PAR.g./(Visc.*exp(PAR.c4.*SnoTemp+PAR.c5.*MR.RHO(:,tn-1)))+PAR.c1.*exp(-PAR.c2.*SnoTemp-PAR.c3.*max(0,MR.RHO(:,tn-1)-RhoNew)))./DQMultCur,RhoMax);
            MR.HS(:,tn) = MR.HS(:,tn-1)./(MR.RHO(:,tn)./MR.RHO(:,tn-1));
          end
          MR.DIA(3,tn) = HS(tn)-sum(MR.HS(:,tn)); % note difference between HSmeas - HSmod after assimilation
          MR.DIA(4,tn) = -DQMultCur; % note assimilation steps required to match HSmeas (negative for decreasing densification)
          %step 3.3.2 if still HSmeas > HSmod (because of RhoMax or because of MAXITER) don't do anything
          if HS(tn) > sum(MR.HS(:,tn))
            % don't do anything
            % later eventually add snow layer MR.DIA(5,tn)
            MR.DIA(5,tn) = -1; % note when densification was too high to meet HS
          end
          %step 3.4 if HSmeas < HSmod reapply densification with gradually increasing densification rate until HSmeas >= HSmod or MR.RHO == RHOmax for all layers
        elseif HS(tn) < sum(MR.HS(:,tn))
          MR.DIA(2,tn) = HS(tn)-sum(MR.HS(:,tn)); % note difference between HSmeas - HSmod before assimilation
          %step 3.4.1 increase deinsification rate
          DQMultCur = 1;
          while mean(MR.RHO(:,tn)) < RhoMax && HS(tn) < sum(MR.HS(:,tn)) && DQMultCur < DQMultMax
            DQMultCur = DQMultCur + DQMultInc;
            MR.RHO(:,tn) = min(MR.RHO(:,tn-1)+MR.RHO(:,tn-1).*PAR.dt.*(MR.OVB(:,tn-1).*PAR.g./(Visc.*exp(PAR.c4.*SnoTemp+PAR.c5.*MR.RHO(:,tn-1)))+PAR.c1.*exp(-PAR.c2.*SnoTemp-PAR.c3.*max(0,MR.RHO(:,tn-1)-RhoNew))).*DQMultCur,RhoMax);
            MR.HS(:,tn) = MR.HS(:,tn-1)./(MR.RHO(:,tn)./MR.RHO(:,tn-1));
          end
          MR.DIA(3,tn) = HS(tn)-sum(MR.HS(:,tn)); % note difference between HSmeas - HSmod after assimilation
          MR.DIA(4,tn) = DQMultCur; % note assimilation steps required to match HSmeas (positive for increasing densification)
          %step 3.4.2 if still HSmeas < HSmod (because of RhoMax or MAXITER) start melting from above
          if HS(tn) < sum(MR.HS(:,tn))
            for lix = size(MR.HS,1):-1:1
              MR.HS(lix,tn) = HS(tn) - sum(MR.HS(1:lix-1,tn));
              if MR.HS(lix,tn) >= 0
                break
              else
                MR.HS(lix,tn) = 0;
              end
            end
            MR.DIA(5,tn) = 1; % note when melt conditions are met
          end
        else %this case should not happen
          error(' ');
        end
        %step 4 recalculate overburden
        nlix = size(MR.HS,1);
        MR.OVB(nlix,tn) = 0;
        for nlix = size(MR.HS,1)-1:-1:1
          MR.OVB(nlix,tn) = sum(MR.HS(nlix+1:end,tn).*MR.RHO(nlix+1:end,tn)./100); % in mm = kg/m2
        end
        %end of loop over days (tn)
      end
      SWE = sum(MR.HS.*MR.RHO./100,1);
      odata(psdix(sepix(tsrix)+1):psdix(sepix(tsrix+1)),serix) = SWE(2:end)';
    end
  end
end
odata(idata==0) = 0;

end