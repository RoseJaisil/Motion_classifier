%data generating file has some dependence 
%change the file
clc,clear
confinded=[20 1200;0 0;0 0]; %large value should be equal to totalTime (below)
free_diff=[0 0;20 1200;0 0];
drift=[0 0;0 0;20 1200];

for i=1
    numParticles=50;%randi([100,200]);
    temp=250;%randi([225,275]);
    volumeEdges=[temp,temp];
    totalTime=1200;%randi([500,700]);
    timeStep=1;
    diffCoefRange=[10,30]; %[randi([8,15]),randi([25,40])];
    confRadRange=[3.5,4];%[randi(3),randi(5)];
    driftVelRange=[0.5*rand,3*rand];
    durationRange=free_diff; %%diffusion, free diffusion, drift
    strictSwitch=0;
    [xCoordMat,yCoordMat,trajClass,errFlag] = simMultiMotionTypeTrajCVMI(numParticles,volumeEdges,totalTime,timeStep,diffCoefRange,confRadRange,driftVelRange,durationRange,strictSwitch);
   % dataVis(xCoordMat,yCoordMat)
end

%# calculate mean standard displacement for all deltaT's
k=1;
%calculating mean square displacement analysis
for ii=1:size(xCoordMat,1);
    data=[xCoordMat(ii,:);yCoordMat(ii,:)];
    data=data';
    nData = size(data,1); %# number of data points
    numberOfdeltaT = floor(nData/4); %in MSD, ii should be 1/4 of number of data points
    msd{ii} = zeros(numberOfdeltaT,3); 
    for j = 1:numberOfdeltaT
        deltaCoords = data(1+j:end,1:2) - data(1:end-j,1:2);
        squaredDisplacement = sum(deltaCoords.^2,2);
        msd{ii}(j,1) = mean(squaredDisplacement); %average
        msd{ii}(j,2) = std(squaredDisplacement); %std
        msd{ii}(j,3) = length(squaredDisplacement); %n
    end
   
end

%mean displacement for all particles at each time step, time average 
for k=1:length(msd);
    j=1:length(msd{1,1});
     mean_values(j,k)=(msd{1,k}(j));
      for m=1:length(mean_values)
          msd_mean(m)=mean(mean_values(m,:));
      end
end
msd_mean=msd_mean';
time=1:timeStep:m;
time=time';

%comparing to different motions
figure,plot(msd_mean);
hold on
title('Free diffusion Fitting Curve');
[free_diff_fit gof_free]=fit(time,msd_mean,'poly1'); %gof=goodness of fit
plot(free_diff_fit)
hold off
figure,plot(msd_mean);
hold on
title('Drift Diffusion Fitting Curve');
[drift_fit gof_drift]=fit(time,msd_mean,'poly2');
plot(drift_fit)
hold off
% gof_free=(gof_free.rsquare);
% gof_drift=(gof_drift.rsquare);
%dataset=[msd_mean;drift_fit];

%slope=linearfit/4;

%classificationLearner
%train classififer with equations 
%use trained classifier for the predition
%load data
%support vector machines
%what is the most ayccurate you choose and go with that
%shouldn't use your brain to classify 
%do something before input into classificationLearner for multi-trajectory
%ground truth (training data) - three types of motion switching

% if gof_free>0.99 || gof_free>gof_drift
%     fprintf('Motion is Free Diffusion')
% elseif gof_drift>0.99 
%     fprintf('Motion is Drift')
% else
%     fprintf('Motion is confined');
% end