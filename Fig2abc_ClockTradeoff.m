% Recovery from jet-lag
clear all;
close all;

% System Parameters
p.w0 = 2*pi; % units of time => 1 day
p.epochT = 45; % for now, let's say this is an integer. This is the number of days that are clean or dirty
p.cleantime = 100; % time before jetlag
p.timedifference = (12/24)*2*pi; % the jump in phase
p.nBug = 6000; % Number of bugs

% Gain parameters
p.kgly = 1/10; % the time scale of gamma going back to rest
p.kmismatch = 3; 
p.gamrest = .3;

% External Noise Statistics
p.sigmaext = .5;
% Internal noise
p.sigmaint = 0; 

% Final params
p.dt = .01;
p.totalT = p.epochT + p.cleantime;

gamlow = .34; 
gamhigh = .73;  


colhighgam = [27 158 119]/256;
collowgam = [117 112 179]/256;  
colk = [217 95 2]/256;


%%%%%%%%%%%%%%%%%%%% High-Low-Kalman - Fig 2b


p.gamma = .34;
[tvec_lag thetamat_lag gammat_lag] = evolveKalmanJetLag(p);
p.gamma = gamlow;
[tvec_low thetamat_low] = evolveBugJetLag(p);
p.gamma = gamhigh;
[tvec_high thetamat_high] = evolveBugJetLag(p);




f=figure;
hold on
plot(tvec_high-p.cleantime,conv(mismatchTime(thetamat_high),ones(1,10)/10,'same'), '--' ,'LineWidth',3,'Color',colhighgam);
plot(tvec_low-p.cleantime,conv(mismatchTime(thetamat_low),ones(1,10)/10,'same'), '--' ,'LineWidth',3,'Color',collowgam);
plot(tvec_lag-p.cleantime,conv(mismatchTime(thetamat_lag),ones(1,10)/10,'same'),'LineWidth',4,'Color',colk);
hold off;
set(gca,'LineWidth',2)
set(gca,'FontSize',12)
set(gca,'FontWeight','bold')
pbaspect([2.5 1 1])
set(gca,'YScale','log')
axis([-10 30 0.3 9])
set(gca,'LineWidth',2)
set(gca,'FontWeight','bold')
f.Renderer = 'Painters';
box on
% saveas(gcf,'Jan19_Error.svg')

f=figure;
plot(tvec_lag-p.cleantime,conv(mean(gammat_lag,2),ones(1,10)/10,'same'),'LineWidth',4,'Color',colk);
hold on
plot(tvec_high-p.cleantime,gamhigh*ones(size(tvec_high)), '--' ,'LineWidth',3,'Color',colhighgam);
plot(tvec_low-p.cleantime,gamlow*ones(size(tvec_low)), '--','LineWidth',3,'Color',collowgam);
hold off;
set(gca,'LineWidth',2)
set(gca,'FontSize',12)
set(gca,'FontWeight','bold')
axis([-10 30 .85*gamlow 1.1*gamhigh])
pbaspect([2.5 1 1])
set(gca,'LineWidth',2)
set(gca,'FontWeight','bold')
f.Renderer = 'Painters';
box on



%%%%%%% Scan $\gamma$, make the tradeoff plot -- 2c

% Compute the time taken and the resting error for the Kalman strategy
p.gamma = .34;
[tvec thetamat gammat] = evolveKalmanRecovery(p);
erroroft = 1-phaseCoherence(thetamat);
varoft = mismatchTime(thetamat);
finalvar = mean(varoft(end-40:end));
finalerror = mean(erroroft(end-40:end)); %average last few points;
minval = finalerror + 0.1*(erroroft(1) - finalerror); % 90% recovery
timetakenKalman = tvec(min(find(conv(1-phaseCoherence(thetamat),ones(1,40)/40,'same') < minval)));
restingKalman = finalvar;
erroroftKalman = erroroft;

p.gamma = gamlow;
[tvec thetamat] = evolveBugRecovery(p);
erroroft = 1-phaseCoherence(thetamat);
varoft = mismatchTime(thetamat);
finalvar = mean(varoft(end-40:end));
finalerror = mean(erroroft(end-40:end)); %average last few points;
minval = finalerror + 0.1*(erroroft(1) - finalerror); % 90% recovery
timetakenLow = tvec(min(find(conv(erroroft,ones(1,40)/40,'same') < minval)));
restingLow = finalvar;


p.gamma = gamhigh;
[tvec thetamat] = evolveBugRecovery(p);
erroroft = 1-phaseCoherence(thetamat);
varoft = mismatchTime(thetamat);
finalvar = mean(varoft(end-40:end));
finalerror = mean(erroroft(end-40:end)); %average last few points;
minval = finalerror + 0.1*(erroroft(1) - finalerror); % 90% recovery
timetakenHigh = tvec(min(find(conv(erroroft,ones(1,40)/40,'same') < minval)));
restingHigh = finalvar;


gammavec = [.28:.005:.95];
timeTaken = [];
restingError = [];
for g = 1:length(gammavec)
	p.gamma = gammavec(g);
	[tvecg thetag] = evolveBugRecovery(p);
	erroroft = 1-phaseCoherence(thetag);
	varoft = mismatchTime(thetag);
	finalvar = mean(varoft(end-40:end));
	finalerror = mean(erroroft(end-40:end)); %average last few points;
	minval = finalerror + 0.1*(erroroft(1) - finalerror); % 90% recovery
	timeTaken(g) = tvec(min(find(conv(1-phaseCoherence(thetag),ones(1,40)/40,'same') < minval)));
	restingError(g) = finalvar;
	disp(howfar({ [g length(gammavec)] }))
end

% % save('tradeoff.mat')

figure;
hold on;
plot(timeTaken,restingError,'o','MarkerSize',3,'LineWidth',1,'Color',[0 0 0]/255,'MarkerFaceColor',[0 0 0]/255)
plot(timetakenKalman,restingKalman,'x','MarkerSize',20,'LineWidth',8,'Color',colk)
plot(timetakenLow,restingLow,'x','MarkerSize',20,'LineWidth',8,'Color',collowgam)
plot(timetakenHigh,restingHigh,'x','MarkerSize',20,'LineWidth',8,'Color',colhighgam)
axis([4 20 0.3 1])
makePretty;
set(gca,'FontSize',20)
set(gca,'FontWeight','normal')
axis([4 20 .35 .7])


%%%%%%% Make a Representative Plot - 2a

p.gamma = 1.5;
p.epochT = 100;
p.nBug = 50;

[tvec_high thetamat_high] = evolveBugJetLag(p);
objtime = @(t) t+ (p.timedifference/(2*pi))*heaviside(t-p.cleantime);

colother = [102,166,30]/256;
sampsize = 25;
f=figure;
for i = 1:sampsize
	plot(tvec_high-p.cleantime,sin(thetamat_high(:,i)),'LineWidth',.5,'Color',colother);
	hold on;
end
hold off
set(gca,'LineWidth',2)
set(gca,'FontWeight','bold')
axis([-4 5 -1.5 1.5])
pbaspect([3 1 1])
f.Renderer = 'Painters';







%%%%% Commonly done tasks

function [popcoh popvar] = calibrationCurve(p, intvec, gammavec)
	popcoh = zeros(length(intvec), length(gammavec));
	popvar = zeros(length(intvec), length(gammavec));
	for inty = 1:length(intvec)
		p.sigmaint = intvec(inty);
		for i = 1:length(gammavec)
			p.gamma = gammavec(i);
			[tvec thetamat] = evolveBugClean(p);
			% keyboard
			popcoht = phaseCoherence(thetamat);
			popvart = var(thetamat,0,2);
			l = length(popvart);
			popvar(inty,i) = mean(popvart(ceil(l/2):end));
			popcoh(inty,i) = mean(popcoht(ceil(l/2):end));
			disp(howfar( { [inty length(intvec)] , [i length(gammavec)] } ))
		end
	end
end


%%%%%%%%%%% Generic Functions


function makePretty()
	set(gca,'LineWidth',2)
	set(gca,'FontSize',20,'FontWeight','normal')
	% set(gca,'FontWeight','bold')
end

% What's the signal at each time?

function siggy = externalSignal(p,t,ps)
	siggy = (sin(2*pi*t)) + (p.sigmaext/sqrt(p.dt))*normrnd(0,1,1,p.nBug);
end

function exty = generateSampleSignal(p,tvec,ps)
	exty = [];
	p.nBug = 1;
	p.dt = tvec(2) - tvec(1); %since this is for visualisation purposes, need to redefine dt as the time step of observation
	for t = 1:length(tvec)
		exty = [exty externalSignal(p,tvec(t),ps)];
	end
end

function cohere = phaseCoherence(thetamat)
	costhet = mean(cos(thetamat),2);
	sinthet = mean(sin(thetamat),2);
	cohere = sqrt(costhet.^2 + sinthet.^2);
end

function cohere = mismatchTime(thetamat)
	costhet = cos(thetamat);
	sinthet = sin(thetamat);

	cosmeanthet = repmat(mean(costhet,2),[1 size(thetamat,2)]);
	sinmeanthet = repmat(mean(sinthet,2),[1 size(thetamat,2)]);

	cosmeanthet = cosmeanthet./(cosmeanthet.^2 + sinmeanthet.^2);
	sinmeanthet = sinmeanthet./(cosmeanthet.^2 + sinmeanthet.^2);


	thetdiff = acos(costhet.*cosmeanthet + sinthet.*sinmeanthet);

	cohere = mean(thetdiff,2)*24/(2*pi);
end

function prc = prcShape(theta)
	prc = cos(theta);
end

function fgam = gamDynamics(p,extsig,theta,gam,t)
	fgam = p.kgly*(-(gam-p.gamrest) - p.kmismatch*(extsig.*sin(theta) - 0.5));
	% fgam = p.kgly*(-(gam-p.gamrest) - p.kmismatch*(extsig-sin(theta)).^2);
end

function stringy = howfar(pairs)
	stringy = '';
	for i = 1:length(pairs)
		stringy = [stringy num2str(pairs{i}(1)) ' of ' num2str(pairs{i}(2)) ' '];
	end
end


%%%%%%%%%% Bugs

%% Basic Time stepping

function [t theta gam]=advanceKalman(p,t,theta,gam,ps)
	dt = p.dt;

	% Advance each bug
	extsig = externalSignal(p,t,ps);
	ftheta = p.w0 + gam.* prcShape(theta) .*( extsig ) + (p.sigmaint/sqrt(dt))*normrnd(0,1,1,p.nBug);
	% fgly = glyDynamicsNoReg(p,extsig,theta,gly,t);
	fgam = gamDynamics(p,extsig,theta,gam,t);
	% keyboard

	theta = theta + dt*ftheta;
	gam = gam + dt*fgam;
	t = t + dt;
end

function [t theta]=advanceBug(p,t,theta,ps)
	dt = p.dt;

	% Advance each bug
	extsig = externalSignal(p,t,ps);
	ftheta = p.w0 + p.gamma.* prcShape(theta) .*( extsig ) + (p.sigmaint/sqrt(dt))*normrnd(0,1,1,p.nBug);
	
	theta = theta + dt*ftheta;
	t = t + dt;
end


%% Protocols


function [tvec thetamat glycmat] = evolveKalmanRecovery(p)
	dt = p.dt;

	innersteps = ceil(.05/dt);% the inner loop -- how many time steps to take before recording the value -- for every n days, do n/dt. needs to be an integer
	outersteps = p.epochT/(innersteps*dt); % The total time of the simulation
	
	% Initialise 
	thetamat = zeros(outersteps + 1,p.nBug); %where i'm going to record things
	theta = 2*pi*(rand(1,p.nBug));
	thetamat(1,:) = theta;

	glycmat = zeros(outersteps + 1,p.nBug); %where i'm going to record things
	
	t = 0;
	tvec = t;
	k = 2;
	gam = p.gamma*ones(size(theta));
	glycmat(1,:) = gam;
	for outer = 1:outersteps
		for inner = 1:innersteps
			% Advance each bug
			[t theta gam]=advanceKalman(p,t,theta,gam,0);
		end
		% Record values
		thetamat(k,:) = theta;
		glycmat(k,:) = gam;
		tvec = [tvec t];
		k = k+1;
	end
end

function [tvec thetamat] = evolveBugRecovery(p)
	dt = p.dt;

	innersteps = ceil(.05/dt);% the inner loop -- how many time steps to take before recording the value -- for every n days, do n/dt. needs to be an integer
	outersteps = p.epochT/(innersteps*dt); % The total time of the simulation
	
	% Initialise 
	thetamat = zeros(outersteps + 1,p.nBug); %where i'm going to record things
	theta = 2*pi*(rand(1,p.nBug));
	thetamat(1,:) = theta;
	
	t = 0;
	tvec = t;
	k = 2;
	for outer = 1:outersteps
		for inner = 1:innersteps
			% Advance each bug
			[t theta]=advanceBug(p,t,theta,0);
		end
		% Record values
		thetamat(k,:) = theta;
		tvec = [tvec t];
		k = k+1;
	end
end


% Jet lag protocols -- start with clean, induce jetlag
function [tvec thetamat glycmat] = evolveKalmanJetLag(p)
	cleantime = ceil(p.cleantime); % how long before jetlag -- needs to be an integer
	dt = p.dt;

	innersteps = ceil(.05/dt);% the inner loop -- how many time steps to take before recording the value -- for every n days, do n/dt. needs to be an integer
	outerstepsRecovery = p.epochT/(innersteps*dt); % The total time of the simulation
	outerstepsClean = cleantime/(innersteps*dt);

	% Initialise 
	thetamat = zeros(outerstepsRecovery + outerstepsClean + 1,p.nBug); %where i'm going to record things
	theta = 0*(rand(1,p.nBug));
	thetamat(1,:) = theta;

	glycmat = zeros(outerstepsRecovery + outerstepsClean + 1,p.nBug); %where i'm going to record things
	
	t = 0;
	tvec = t;
	k = 2;
	gam = p.gamma*ones(size(theta));
	glycmat(1,:) = gam;
	% Clean period
	for outer = 1:outerstepsClean
		for inner = 1:innersteps
			% Advance each bug
			[t theta gam]=advanceKalman(p,t,theta,gam,0);
		end
		% Record values
		thetamat(k,:) = theta;
		glycmat(k,:) = gam;
		tvec = [tvec t];
		k = k+1;
	end

	% Jet lag
	theta = theta + p.timedifference*rand(1,p.nBug).*((randi([0 1],1,p.nBug) - 0.5)*2);

	% Recovery period
	for outer = 1:outerstepsRecovery
		for inner = 1:innersteps
			% Advance each bug
			[t theta gam]=advanceKalman(p,t,theta,gam,0);
		end
		% Record values
		thetamat(k,:) = theta;
		glycmat(k,:) = gam;
		tvec = [tvec t];
		k = k+1;
	end
end


function [tvec thetamat] = evolveBugJetLag(p)
	cleantime = ceil(p.cleantime); % how long before jetlag -- needs to be an integer
	dt = p.dt;

	innersteps = ceil(.05/dt);% the inner loop -- how many time steps to take before recording the value -- for every n days, do n/dt. needs to be an integer
	outerstepsRecovery = p.epochT/(innersteps*dt); % The total time of the simulation
	outerstepsClean = cleantime/(innersteps*dt);

	% Initialise 
	thetamat = zeros(outerstepsRecovery + outerstepsClean + 1,p.nBug); %where i'm going to record things
	theta = 0.1*(rand(1,p.nBug)) + 0.2;
	thetamat(1,:) = theta;
	
	t = 0;
	tvec = t;
	k = 2;
	% Clean period
	for outer = 1:outerstepsClean
		for inner = 1:innersteps
			% Advance each bug
			[t theta]=advanceBug(p,t,theta,0);
		end
		% Record values
		thetamat(k,:) = theta;
		tvec = [tvec t];
		k = k+1;
	end

	% Jet lag
	theta = theta + 2*pi*p.timedifference*rand(1,p.nBug).*((randi([0 1]) - 0.5)*2);

	% Recovery period
	for outer = 1:outerstepsRecovery
		for inner = 1:innersteps
			% Advance each bug
			[t theta]=advanceBug(p,t,theta,0);
		end
		% Record values
		thetamat(k,:) = theta;
		tvec = [tvec t];
		k = k+1;
	end
end




function [tvec thetamat glycmat] = evolveKalmanClean(p)
	dt = p.dt;

	innersteps = ceil(.05/dt);% the inner loop -- how many time steps to take before recording the value -- for every n days, do n/dt. needs to be an integer
	outersteps = p.epochT/(innersteps*dt); % The total time of the simulation

	% Initialise 
	thetamat = zeros(p.numEp*2*outersteps + 1,p.nBug); %where i'm going to record things
	theta = (rand(1,p.nBug));
	thetamat(1,:) = theta;

	glycmat = zeros(p.numEp*2*outersteps + 1,p.nBug); %where i'm going to record things
	
	t = 0;
	tvec = t;
	k = 2;
	gam = p.gamma*ones(size(theta));
	glycmat(1,:) = gam;
	for ep=1:2*p.numEp
		% First, the clean epoch
		for outer = 1:outersteps
			for inner = 1:innersteps
				% Advance each bug
				[t theta gam]=advanceKalman(p,t,theta,gam);
			end
			% Record values
			thetamat(k,:) = theta;
			glycmat(k,:) = gam;
			tvec = [tvec t];
			k = k+1;
		end
	end
end

function [tvec thetamat glycmat] = evolveKalmanDirty(p)
	dt = p.dt;

	innersteps = ceil(.05/dt);% the inner loop -- how many time steps to take before recording the value -- for every n days, do n/dt. needs to be an integer
	outersteps = p.epochT/(innersteps*dt); % The total time of the simulation

	% dirty and clean
	cleanint = p.sigmaint;
	dirtyint = p.sigmaintDirty;
	cleanext = p.sigmaext;
	dirtyext = p.sigmaextDirty;

	% Initialise 
	thetamat = zeros(p.numEp*2*outersteps + 1,p.nBug); %where i'm going to record things
	theta = (rand(1,p.nBug));
	thetamat(1,:) = theta;

	glycmat = zeros(p.numEp*2*outersteps + 1,p.nBug); %where i'm going to record things
	

	t = 0;
	tvec = t;
	k = 2;
	gam = p.gamma*ones(size(theta));
	glycmat(1,:) = gam;
	for ep=1:p.numEp
		% First, the clean epoch
		p.sigmaint = cleanint;
		p.sigmaext = cleanext;
		for outer = 1:outersteps
			for inner = 1:innersteps
				% Advance each bug
				[t theta gam]=advanceKalman(p,t,theta,gam);
			end
			% Record values
			thetamat(k,:) = theta;
			glycmat(k,:) = gam;
			tvec = [tvec t];
			k = k+1;
		end

		% Then, the dirty epoch
		p.sigmaint = dirtyint;
		p.sigmaext = dirtyext;
		for outer = 1:outersteps
			for inner = 1:innersteps
				% Advance each bug
				[t theta gam]=advanceKalman(p,t,theta,gam);
				theta = theta + p.rotationSize*rand(1,p.nBug).*((rand(1,p.nBug) < p.rotRate*dt).*(randi([0 1]) - 0.5)*2);
			end
			% Record values
			thetamat(k,:) = theta;
			glycmat(k,:) = gam;
			tvec = [tvec t];
			k = k+1;
		end
	end
end


%%%%% Non-Kalman

function [tvec thetamat] = evolveBugClean(p)
	dt = p.dt;

	innersteps = ceil(.05/dt);% the inner loop -- how many time steps to take before recording the value -- for every n days, do n/dt. needs to be an integer
	outersteps = p.epochT/(innersteps*dt); % The total time of the simulation

	% Initialise 
	thetamat = zeros(p.numEp*2*outersteps + 1,p.nBug); %where i'm going to record things
	theta = 2*pi*(0.2*rand(1,p.nBug) - 0.1);
	thetamat(1,:) = theta;

	t = 0;
	tvec = t;
	k = 2;
	for ep=1:2*p.numEp
		% First, the clean epoch
		for outer = 1:outersteps
			for inner = 1:innersteps
				% Advance each bug
				[t theta]=advanceBug(p,t,theta);
			end
			% Record values
			thetamat(k,:) = theta;
			tvec = [tvec t];
			k = k+1;
		end
	end
end

function [tvec thetamat] = evolveBugDirty(p)
	% Evolve a set of bugs over clean and dirty environments. For 
	dt = p.dt;

	innersteps = ceil(.1/dt);% the inner loop -- how many time steps to take before recording the value -- for every n days, do n/dt. needs to be an integer
	outersteps = p.epochT/(innersteps*dt); % The total time of the simulation

	% dirty and clean
	cleanint = p.sigmaint;
	dirtyint = p.sigmaintDirty;
	cleanext = p.sigmaext;
	dirtyext = p.sigmaextDirty;


	% Initialise 
	thetamat = zeros(p.numEp*2*outersteps + 1,p.nBug); %where i'm going to record things
	theta = 2*pi*(0.2*rand(1,p.nBug) - 0.1);
	thetamat(1,:) = theta;
	t = 0;
	tvec = t;
	k = 2;
	for ep = 1:p.numEp
		% First, the clean epoch
		p.sigmaint = cleanint;
		p.sigmaext = cleanext;
		for outer = 1:outersteps
			for inner = 1:innersteps
				% Advance each bug
				[t theta]=advanceBug(p,t,theta);
			end
			% Record values
			thetamat(k,:) = theta;
			tvec = [tvec t];
			k = k+1;
		end

		% Now, the noisy epoch
		p.sigmaint = dirtyint;
		p.sigmaext = dirtyext;
		for outer = 1:outersteps
			for inner = 1:innersteps
				% Advance each bug
				[t theta]=advanceBug(p,t,theta);
				theta = theta + p.rotationSize*rand(1,p.nBug).*((rand(1,p.nBug) < p.rotRate*dt).*(randi([0 1]) - 0.5)*2);
			end
			% Record values
			thetamat(k,:) = theta;
			tvec = [tvec t];
			k = k+1;
		end
	end
end














