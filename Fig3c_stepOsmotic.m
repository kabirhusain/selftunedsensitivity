clear all;
close all;

% System Parameters
p.tau_f = 3; % set units of time
p.tau_s = 6; % the slow responder
p.tau = 1; % fixed timescale
p.adaptationTimeScale = 6;
p.kappaThresh = 1;
p.y0 = 0;

p.P0 = 10;
p.P1 = 20;
p.sigmaext = 1;

meansig = @(tvec) p.P0*double(tvec <= 0) + p.P1*double(tvec > 0);

% Simulation Parameters
p.totalT = 50;
p.dt = .01;
p.nBug = 6000;

colhighgam = [27 158 119]/256;
collowgam = [117 112 179]/256;
colk = [217 95 2]/256;

%%%%%% The slow, fast, and adaptive strategies.
p.tau = p.tau_f;
[tvechigh glymathigh] = evolvePressureJumpSimple(p);
p.tau = p.tau_s;
[tveclow glymatlow] = evolvePressureJumpSimple(p);
[tvec glymat kappamat alphmat] = evolvePressureJumpKalman(p);

[tRecoveryK restingErrorK]=ParetoPoints(p,tvec,glymat);
[tRecovery_high restingError_high]=ParetoPoints(p,tvechigh,glymathigh);
[tRecovery_low restingError_low]=ParetoPoints(p,tveclow,glymatlow);

%%%%%%%%%% Plot examples of response to a jump in pressure

% plot osmotic pressure
figure;
plot(tvec,meansig(tvec),'-k','LineWidth',3)
axis([-5 22 .8*p.P0 1.1*(p.P1)])
makePretty
pbaspect([4 1 1])
box on;

% plot pressure mismatch
figure;
meany = meansig(tvec)' - mean(glymat,2);
stdy = std((meansig(tvec)' - glymat),[],2);
hold on;
plot(tvec, meansig(tvec)' - mean(glymathigh,2),'--','LineWidth',3,'Color',colhighgam)
plot(tvec, meansig(tvec)' - mean(glymatlow,2),'--','LineWidth',3,'Color',collowgam)
plot(tvec, meansig(tvec)' - mean(glymat,2),'-','LineWidth',3,'Color',colk)
hold off
axis([-5 22 -2 12])
pbaspect([2 1 1])
makePretty
box on

figure;
hold on;
plot(tvec, mean(kappamat,2),'-','LineWidth',3,'Color',colk)
plot(tvec, ones(size(tvec)),'--','LineWidth',3,'Color',colhighgam)
plot(tvec, zeros(size(tvec)),'--','LineWidth',3,'Color',collowgam)
hold off
axis([-5 22 -.2 1.2])
pbaspect([3 1 1])
makePretty
box on


%%%%%%%%%%%%%%%%%%%%%%% Functions

function shadyPlot(tvy,meany,stdy,col);
	tvy = tvy(:)';
	meany = meany(:)';
	stdy = stdy(:)';

	hold on;

	% keyboard

	efy=fill([tvy fliplr(tvy)],[meany+stdy fliplr(meany-stdy)],col);
	alpha(efy,.5)
	% set(efy,'LineWidth',0)
	stairs(tvy,meany,'Color',col,'LineWidth',3)
end


function shadyPlotD(tvy,meany,stdy,col);
	tvy = tvy(:)';
	meany = meany(:)';
	stdy = stdy(:)';

	hold on;

	% keyboard

	efy=fill([tvy fliplr(tvy)],[meany+stdy fliplr(meany-stdy)],col);
	alpha(efy,.5)
	% set(efy,'LineWidth',0)
	stairs(tvy,meany,'--','Color',col,'LineWidth',3)
end



%%%%%%%% Compute Resting Error
function [tRecovery restingError]=ParetoPoints(p,tvec,glymat)
	tRecovery = tvec(min(intersect(find((p.P1 - mean(glymat,2))/(p.P1-p.P0) < .2),find(tvec>0))));
	tosample = max(find(tvec<-0.5));
	stdGly = std(glymat,0,2);
	restingError = mean(stdGly(tosample-20:tosample));
end

%%%%%%%% Euler Steps

function [t gly kappa alph] = advanceYeastKalman(p,t,gly,alph,extsig)
	dt = p.dt;
	c = 2/p.adaptationTimeScale;

	% kappa = double(alph-p.y0 > p.kappaThresh);
	kappa = ((alph).^4)./((p.kappaThresh+p.y0).^4+(alph).^4);

	fgly = -(1 - kappa).*(gly - extsig)/p.tau_s - kappa.*(gly-extsig)/p.tau_f;
	falpha = -(alph - (p.y0 + extsig - gly))/p.adaptationTimeScale;
	% fx = (adaptory-p.y0);
	% fy = -adaptorx*(c^2)/4 - c*(adaptory-p.y0) + (extsig-gly);

	% adaptorx = adaptorx + dt*fx;
	% adaptory = adaptory + dt*fy;
	alph = alph + dt*falpha;
	gly = gly + dt*fgly;
	t = t + dt;
end

function [t gly] = advanceYeastSimple(p,t,gly,extsig)
	dt = p.dt;

	fgly = -(gly - extsig)/p.tau;

	gly = gly + dt*fgly;
	t = t + dt;
end



%%%%%%%%%%%%%%%%%%% Protocols

function [tvec glymat] = evolvePressureJumpSimple(p)
	dt = p.dt;

	innersteps = ceil(.05/dt);% the inner loop -- how many time steps to take before recording the value -- for every n days, do n/dt. needs to be an integer
	outersteps = 2*p.totalT/(innersteps*dt); % The total time of the simulation
	
	% Initialise 
	glymat = zeros(outersteps + 1,p.nBug); %where i'm going to record things
	gly = p.P0 + 0*(rand(1,p.nBug));
	glymat(1,:) = gly;

	t = -p.totalT;
	tvec = t;
	k=2;
	for outer = 1:outersteps
		for inner = 1:innersteps
			extsig = p.P0*double(t < 0) + p.P1*double(t > 0) + (p.sigmaext/sqrt(dt))*normrnd(0,1,1,p.nBug);
			% keyboard;
			[t gly] = advanceYeastSimple(p,t,gly,extsig);
		end
		glymat(k,:) = gly;
		tvec = [tvec t];
		k = k+1;
	end
end

function [tvec glymat kappamat alphmat] = evolvePressureJumpKalman(p)
	dt = p.dt;

	innersteps = ceil(.05/dt);% the inner loop -- how many time steps to take before recording the value -- for every n days, do n/dt. needs to be an integer
	outersteps = 2*p.totalT/(innersteps*dt); % The total time of the simulation
	
	% Initialise 
	glymat = zeros(outersteps + 1,p.nBug); %where i'm going to record things
	gly = p.P0 + 0*(rand(1,p.nBug));
	glymat(1,:) = gly;

	kappamat = zeros(outersteps + 1,p.nBug); %where i'm going to record things
	kappa = zeros(1,p.nBug);
	kappamat(1,:) = kappa;

	alphmat = zeros(outersteps + 1,p.nBug); %where i'm going to record things
	alph = zeros(1,p.nBug);
	alphmat(1,:) = alph;

	t = -p.totalT;
	tvec = t;
	k=2;
	for outer = 1:outersteps
		for inner = 1:innersteps
			extsig = p.P0*double(t <= 0) + p.P1*double(t > 0) + (p.sigmaext/sqrt(dt))*normrnd(0,1,1,p.nBug);
			% keyboard;
			[t gly kappa alph] = advanceYeastKalman(p,t,gly,alph,extsig);
		end
		glymat(k,:) = gly;
		kappamat(k,:) = kappa;
		alphmat(k,:) = alph;
		tvec = [tvec t];
		k = k+1;
	end
end


function [tvec glymat kappamat adaptorxmat adaptorymat] = evolveManyPressureJumpKalman(p)
	dt = p.dt;

	epnum = 4;
	pdiff = p.P1-p.P0;

	innersteps = ceil(.05/dt);% the inner loop -- how many time steps to take before recording the value -- for every n days, do n/dt. needs to be an integer
	outersteps = 2*p.totalT/(innersteps*dt); % The total time of the simulation
	
	% Initialise 
	glymat = zeros(epnum*outersteps + 1,p.nBug); %where i'm going to record things
	gly = p.P0 + 0*(rand(1,p.nBug));
	glymat(1,:) = gly;

	kappamat = zeros(epnum*outersteps + 1,p.nBug); %where i'm going to record things
	kappa = zeros(1,p.nBug);
	kappamat(1,:) = kappa;

	adaptoxmat = zeros(epnum*outersteps + 1,p.nBug); %where i'm going to record things
	adaptorx = zeros(1,p.nBug);
	adaptoxmat(1,:) = adaptorx;

	adaptoymat = zeros(epnum*outersteps + 1,p.nBug); %where i'm going to record things
	adaptory = zeros(1,p.nBug);
	adaptoymat(1,:) = adaptory;	

	t = 0;
	tvec = t;
	k=2;
	for ep = 1:epnum
		pnow = p.P0 + (ep-1)*pdiff;
		for outer = 1:outersteps
			for inner = 1:innersteps
				extsig = pnow + (p.sigmaext/sqrt(dt))*normrnd(0,1,1,p.nBug);
				% keyboard;
				[t gly kappa adaptorx adaptory] = advanceYeastKalman(p,t,gly,adaptorx,adaptory,extsig);
			end
			glymat(k,:) = gly;
			kappamat(k,:) = kappa;
			adaptorxmat(k,:) = adaptorx;
			adaptorymat(k,:) = adaptory;
			tvec = [tvec t];
			k = k+1;
		end
	end
end


%%%%%%%% Generic Functions
function stringy = howfar(pairs)
	stringy = '';
	for i = 1:length(pairs)
		stringy = [stringy num2str(pairs{i}(1)) ' of ' num2str(pairs{i}(2)) ' '];
	end
end

function makePretty()
	set(gca,'LineWidth',2)
	set(gca,'FontSize',12)
	set(gca,'FontWeight','bold')
end