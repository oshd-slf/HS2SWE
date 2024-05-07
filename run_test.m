% Clear

clear; home;

% Load data used in hnw-study

load('D:\slf_code\oi\data\point_swe_hs.mat','hs_stake','swe_hs2swe')

% Run standalone HS2SWE model

swe_test = HS2SWE(hs_stake);

% Plot results

plot(swe_hs2swe(:),swe_test(:),'.')
xlabel("SWE REF")
ylabel("SWE TEST")

% Assert results

assert(max(abs(swe_hs2swe(:)-swe_test(:)))<1e-10,"Mismatch in results")


%% Extract an example time series from Weissfluhjoch

% Clear

clear; home;

% Load data used in hnw-study

load('D:\slf_code\oi\data\point_swe_hs.mat','acro','time','hs_stake','swe_hs2swe','swe_profile')

% Get data for Weisfluhjoch

istat = find(acro=="SLF.5WJ");
itime = find(time >= datenum(2016,9,1) & time <= datenum(2022,9,1));

time = time(itime);
hs_obs = hs_stake(itime,istat);

% Run model and also extract verification data

swe_test = HS2SWE(hs_obs);
swe_ref = swe_hs2swe(itime,istat);

% Plot results

plot(swe_ref(:),swe_test(:),'.')
xlabel("SWE REF")
ylabel("SWE TEST")

% Assert results

assert(max(abs(swe_test(:)-swe_ref(:)))<1e-10,"Mismatch in results")

% Save data to table

time = datestr(time,"yyyy-mm-dd HH:MM");
swe_obs = swe_profile(itime,istat);

data_weissfluhjoch = table(time,hs_obs,swe_obs,swe_ref);

writetable(data_weissfluhjoch)
