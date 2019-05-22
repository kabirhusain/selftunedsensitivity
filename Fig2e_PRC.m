clear all;
close all;

p.w0 = 2*pi; % units of time => 1 day
p.epochT = 30; % for now, let's say this is an integer. This is the number of days that are clean or dirty
p.numEp = 2; % Number of clean-dirty pairs
p.nBug = 1; % Number of bugs -- don't worry, look after
repli = 2000;

% Glycogen parameters
p.kglyMake = 16;
p.kglyConsume = .4;

% Dirty period parameters
p.rotRate = 1; % average rotations per bug per day
p.rotMax = 2*pi*(12/24); % rotation goes between -rotMax and rotMax; sigma_int is 1/2 of this

% Relaxation to the limit cycle: 
p.alpha = 10;
p.dt = .001;
p.gamma = .5;

% PRC variables
p.PRCsteps = 100;
p.PRCpulse = (5/24);

colhighgam = [27 158 119]/256;
collowgam = [117 112 179]/256;  
colk = [217 95 2]/256;

% First entrain the cell so that glycogen gets to its resting `high' level.
[tvec xmat ymat glycmat] = evolveKalmanClean(p);
[t0b deltatheta] = generatePRC(p,xmat(end),ymat(end),glycmat(end));
PRC_before = 24*deltatheta/(2*pi);

% Now do a dirty period

epvec = [.1:.05:15];
sigvec = [0.2:.05:6];

p.nBug = repli; % now we do many bugs = many realisations of dirty periods
for sig = 1:length(sigvec)
	prcmax = zeros(size(epvec));
	for ep = 1:length(epvec)
		epochT = epvec(ep);
		rotMax = sigvec(sig);
		[tvec xmatD ymatD glycmatD] = evolveKalmanDirty(p,epochT,rotMax,repmat(xmat(end),[1,p.nBug]),repmat(ymat(end),[1,p.nBug]),repmat(glycmat(end),[1,p.nBug]));
		[t0at deltathetat] = generatePRC(p,xmatD(end,:),ymatD(end,:),glycmatD(end,:));

		prcmax(ep) = 24*mean(max(deltathetat,[],1))/(2*pi);

		disp(howfar( { [ep length(epvec)], [sig length(sigvec)] } ))
	end
	t0at = [];
	deltathetat = [];
	prcMaxVal{sig} = prcmax;
end

% Convert to sigma_int -- average magnitude of phase kick
sigvec = (sigvec/2)*(24/(2*pi));

prcHeight = zeros(length(sigvec),length(epvec));
for sig = 1:length(sigvec)
	prcHeight(sig,:) = prcMaxVal{sig};
end

figure;
imagesc(epvec,sigvec,prcHeight/max(PRC_before));
set(gca,'YDir','normal')
colormap('jet')
h=colorbar;
makePretty

%%%% Functions
function fgly = glyDynamics(p,extsig,x,y,gly,t)
	% keyboard
	fgly = p.kglyMake*((extsig>0).*(y>0)) - p.kglyConsume*((gly).*(y<0));
end

function gam = gly2gam(gly,p)
	deltaatp = max(0,-1.6*gly + 72);
	gam = (.0123*deltaatp + .2607); %.*(y>0) + p.daytimegamma.*(y<0); %day time coupling is a constant
end

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

function [x y] = rotateTheRotters(x,y,p,gam,extsig,t)
	% Rotate some bugs -- let's avoid another loop
	phi = 2*p.rotMax*(rand(1,p.nBug) - 0.5);
	phi = phi.*(rand(1,p.nBug) < p.rotRate*p.dt);

	realsig = 1-extsig;
	x0 = gam.*realsig;
	x = x-x0;	

	xn = x.*cos(phi) - y.*sin(phi);
	yn = y.*cos(phi) + x.*sin(phi);
	x = xn + x0;
	y = yn;
end

%%%% Time-stepping

function [t x y gly gam]=advanceKalman(p,t,x,y,gly,gam,extsig)
	dt = p.dt;

	realsig = 1 - extsig; % 1 when night, 0 when day

	% Decide x0 as a function of the external signal
	x0 = gam.*(realsig);

	fx = -p.w0*y + p.alpha*(x-x0).*(1 - y.^2 - (x-x0).^2);
	 % + (p.sigmaint/sqrt(dt))*normrnd(0,1,1,p.nBug);
	fy = p.w0*(x-x0) + p.alpha*y.*(1 - (x-x0).^2 - y.^2);
	 % + (p.sigmaint/sqrt(dt))*normrnd(0,1,1,p.nBug);

	fgly = glyDynamics(p,extsig,x,y,gly,t);
	% keyboard

	x = x + dt*fx;
	y = y + dt*fy;
	gly = gly + dt*fgly;
	gam = gly2gam(gly,p);
	t = t + dt;
end


%%%%%%% Protocols

function [tvec xmat ymat glycmat] = evolveKalmanClean(p)
	dt = p.dt;

	innersteps = ceil(.05/dt);% the inner loop -- how many time steps to take before recording the value -- for every n days, do n/dt. needs to be an integer
	outersteps = p.epochT/(innersteps*dt); % The total time of the simulation

	% Initialise 
	xmat = zeros(outersteps + 1,p.nBug); %where i'm going to record things
	ymat = zeros(outersteps + 1,p.nBug); %where i'm going to record things
	
	theta = 2*pi*(rand(1,p.nBug));
	x = cos(theta);
	y = sin(theta);
	
	xmat(1,:) = x;
	ymat(1,:) = y;

	glycmat = zeros(outersteps + 1,p.nBug); %where i'm going to record things
	gly = 15 + .2*(rand(1,p.nBug));
	glycmat(1,:) = gly;

	t = 0;
	tvec = t;
	k = 2;
	gam = p.gamma*ones(size(gly));
	
	for outer = 1:outersteps
		for inner = 1:innersteps
			[t x y gly gam]=advanceKalman(p,t,x,y,gly,gam,sin(2*pi*t)>0);
		end
		% Record values
		xmat(k,:) = x;
		ymat(k,:) = y;
		glycmat(k,:) = gly;
		tvec = [tvec t];
		k = k+1;
	end
end

function [tvec xmat ymat glycmat] = evolveKalmanDirty(p,epochT,rotMax,x,y,gly)
	dt = p.dt;

	p.epochT = epochT;
	p.rotMax = rotMax;

	innersteps = ceil(.05/dt);% the inner loop -- how many time steps to take before recording the value -- for every n days, do n/dt. needs to be an integer
	outersteps = round(p.epochT/(innersteps*dt)); % The total time of the simulation

	% Initialise 
	xmat = zeros(outersteps + 1,p.nBug); %where i'm going to record things
	ymat = zeros(outersteps + 1,p.nBug); %where i'm going to record things
		
	xmat(1,:) = x;
	ymat(1,:) = y;

	glycmat = zeros(outersteps + 1,p.nBug); %where i'm going to record things
	glycmat(1,:) = gly;

	t = 0;
	tvec = t;
	gam = gly2gam(gly);
	k = 2;
	for outer = 1:outersteps
		for inner = 1:innersteps
			[t x y gly gam]=advanceKalman(p,t,x,y,gly,gam,sin(2*pi*t)>0);
			[x y] = rotateTheRotters(x,y,p,gam,sin(2*pi*t)>0,t);
		end
		% Record values
		xmat(k,:) = x;
		ymat(k,:) = y;
		glycmat(k,:) = gly;
		tvec = [tvec t];
		k = k+1;
	end
end


function [theta0 deltatheta] = generatePRC(p,x,y,gly)
	xorig = x;
	yorig = y;
	glyorig = gly;

	% Go through the PRC
	deltatheta = zeros(p.PRCsteps,size(x,2));
	theta0 = zeros(p.PRCsteps,size(x,2));
	xp = [];
	yp = [];
	for pr = 1:p.PRCsteps
		waittime = (pr-1)/p.PRCsteps;
		
		% delay in white light
		t = 0;
		x = xorig;
		y = yorig;
		gly = glyorig;
		gam = gly2gam(gly);
		while t < waittime
			[t x y gly gam]=advanceKalman(p,t,x,y,gly,gam,1); % in light
		end

		% save these as initial conditions for what follows
		xinit = x;
		yinit = y;
		glyinit = gly;

		% Now do five hours in dark followed by 10 hours in day
		t = 0;
		x = xinit;
		y = yinit;
		gly = glyinit;
		gam = gly2gam(gly);
		while t < p.PRCpulse
			[t x y gly gam]=advanceKalman(p,t,x,y,gly,gam,0); % in dark
		end
		t = 0;
		while t < 2*p.PRCpulse
			[t x y gly gam]=advanceKalman(p,t,x,y,gly,gam,1); % in light
		end
		% final values -- pulse ones:
		xpulse = x;
		ypulse = y;

		% Now do 15 in day as a reference
		t = 0;
		x = xinit;
		y = yinit;
		gly = glyinit;
		gam = gly2gam(gly);
		while t < 3*p.PRCpulse
			[t x y gly gam]=advanceKalman(p,t,x,y,gly,gam,1); % in light
		end

		% keyboard

		% Compute the delta theta
		% keyboard
		for rep = 1:size(x,2)
			maggy = sqrt(y(rep)^2 + x(rep)^2); %should be one, but lets just be safe
			thetaf = atan2(y(rep)/maggy,x(rep)/maggy);
			rotmat = [cos(thetaf) sin(thetaf); -sin(thetaf) cos(thetaf)];
			finalv = rotmat*[xpulse(rep); ypulse(rep)];
			finalv = finalv/sqrt(finalv'*finalv);
			deltatheta(pr,rep) = atan2(finalv(2),finalv(1));
			theta0(pr,rep) = thetaf;
		end
	end
	% keyboard
end

