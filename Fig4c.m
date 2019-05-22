clear all;
close all;


colhighgam = [27 158 119]/256;
collowgam = [117 112 179]/256;
colk = [217 95 2]/256;

gamlow = .1;
gamhigh = .2;

p.vel = 0.1;
p.velMult = 1;
p.sigmaext = .5;
p.sigmaint = 0.01;
p.kappa = 0.2;
p.taug = 100;
p.gamma = .08; % fixed gamma

p.sigClean = p.sigmaint;
p.sigDirty = 1;

p.nBug = 600;
p.totalT = 1500; % actually the epoch time
p.cleanTime = 800;
p.dirtyTime = 100;
p.epochNum = 1;
p.dt = 0.05;


%%%%% Plot tracking error 

lowgam = steadyStateGamma(p.sigClean, p.sigmaext, p.kappa)
highgam = steadyStateGamma(p.sigDirty, p.sigmaext, p.kappa)


p.dirtyTime = 500;
p.gamma = lowgam;
low_results = ballFixed(p);
tv = low_results.tvec - p.cleanTime;
indirty = intersect(find(tv>0), find(tv<p.dirtyTime));
trackerror = mean((low_results.xmat - low_results.realxmat).^2,2);
lowgamma_perf = mean(trackerror(indirty));

p.gamma = highgam;
high_results = ballFixed(p);
tv = high_results.tvec - p.cleanTime;
indirty = intersect(find(tv>0), find(tv<p.dirtyTime));
trackerror = mean((high_results.xmat - high_results.realxmat).^2,2);
highgamma_perf = mean(trackerror(indirty));


p.gamma = lowgam;
dirtvec = [100 250 500 1000];

taugvec = 10.^[0.8:.1:5.5];

ditty = 1;
for ditty = 1:length(dirtvec)
	kal_perf{ditty} = [];
	p.dirtyTime = dirtvec(ditty);
	for tt = 1:length(taugvec)
		p.taug = taugvec(tt);
		kr = ballKalman(p);

		kal_results = kr;

		tv = kal_results.tvec - p.cleanTime;
		indirty = intersect(find(tv>0), find(tv<p.dirtyTime));
		trackerror = mean((kal_results.xmat - kal_results.realxmat).^2,2);
		kal_perf{ditty} = [kal_perf{ditty} mean(trackerror(indirty))];

		disp(howfar({[tt length(taugvec)],[ditty length(dirtvec)]}))
	end
end

save('trackerror.mat');


load('trackerror.mat');

figure;
hold on
plot([taugvec(1) taugvec(end)],[highgamma_perf highgamma_perf],'--','LineWidth',4,'Color',colhighgam)
plot([taugvec(1) taugvec(end)],[lowgamma_perf lowgamma_perf],'--','LineWidth',4,'Color',collowgam)
cols = jet(length(dirtvec));
markers = {'x','o','s','+'};
for ditty = 1:length(dirtvec)
	plot(taugvec,kal_perf{ditty},markers{ditty},'MarkerSize',8,'LineWidth',1,'Color',colk)
	plot(taugvec,kal_perf{ditty},'-','MarkerSize',10,'LineWidth',1,'Color',colk)
end
hold off
axis([taugvec(1) taugvec(end) 1.2 6])
set(gca,'XScale','log')
% set(gca,'YScale','log')
box on;
pbaspect([1 1 1])
makePretty
saveas(gcf,'svg_trackingError.svg')

tauvec = 10.^[-3:.1:3];
exacto = analyticPerf((1/rateGamma(p.sigDirty, p.sigmaext, p.kappa))*tauvec,p,lowgam,highgam);

figure;
hold on
plot([0.007 3000],[highgamma_perf highgamma_perf],'--','LineWidth',4,'Color',colhighgam)
plot([0.007 3000],[lowgamma_perf lowgamma_perf],'--','LineWidth',4,'Color',collowgam)
cols = jet(length(dirtvec));
markers = {'x','o','s','+'};
for ditty = 1:length(dirtvec)
	plot(taugvec/dirtvec(ditty),kal_perf{ditty},markers{ditty},'MarkerSize',6,'LineWidth',1,'Color',colk)
end
plot(tauvec,exacto,'-','LineWidth',1.5,'Color',[0 0 0]/255);

hold off
axis([0.0063 1000 1.2 6])
set(gca,'XScale','log')
% set(gca,'YScale','log')
box on;
pbaspect([1.5 1 1])
makePretty
saveas(gcf,'svg_trackingError_scaled.svg')



%%%%%% Functions

function ratey = rateGamma(sigint, sigext, kappa)
	fung = @(g) g^2 - (kappa^2)*(2/pi)*((g^2 +1)*(sigext^2) + (sigint^2)*((g-1)^2)/(1-(g-1)^2));
	gam = fzero(fung,0.2);

	nume = ((gam-2)^2)*(gam^3)*(sigext^2) + (gam-1)*sigint^2;
	deno = ((gam-2)^2)*(gam^2)*sqrt((gam^2 + 1)*(sigext^2) + (sigint^2)*(((1-gam)^2)/(1-((1-gam)^2))) );
	ratey = 1 - kappa*sqrt(2/pi)*nume/deno;
end

function gam = steadyStateGamma(sigint, sigext, kappa)
	fung = @(g) g^2 - (kappa^2)*(2/pi)*((g^2 +1)*(sigext^2) + (sigint^2)*((g-1)^2)/(1-(g-1)^2));
	gam = fzero(fung,0.2);
end

function exacto = analyticPerf(t,p,gamlow,gamhigh)
	g0 = gamlow;
	g1 = gamhigh;
	g0mg1 = g0-g1;
	g1mg0 = g1-g0;
	k = t.^(-1);

	coeffExt = (p.sigmaext^2)*g0mg1*t.*( ((g1^2)/g0mg1)*k + (g0 + 3*g1)/2 - 0.5*exp(-2*k).*(g0 + g1*( 4*exp(k) - 1  ) ) );

	coeffInt = -(p.sigDirty^2)* ( (g1-1)^2 + t.*( g1*( atanh(1 - g1 + g1mg0*exp(-k)) - atanh(1-g0) ) + log( (1-(g1/g0))*exp(-k) + (g1/g0) ) ) ) /(g1*(g1-2));

	exacto = coeffInt + coeffExt;
end

function [x gam realx t] = advanceBallKalman(p,x,gam,realx,t)
	% measured value
	realx = realx + p.vel;
	meas = realx + (p.sigmaext)*normrnd(0,1,1,p.nBug);
	pred = x + p.vel + p.sigmaint*normrnd(0,1,1,p.nBug);
	
	x = pred + gam.*(meas - pred);

	% new gamma
	gam = min(gam*(1 - (1/p.taug)) + (p.kappa/p.taug)*abs(x - meas),1);
	t = t+1;
end

function [x realx t] = advanceBallFixed(p,x,realx,t)
	% measured value
	realx = realx + p.vel;
	meas = realx + (p.sigmaext)*normrnd(0,1,1,p.nBug);
	pred = x + p.vel + p.sigmaint*normrnd(0,1,1,p.nBug);
	
	x = pred + p.gamma*(meas - pred);

	% new gamma
	t = t+1;
end


function simresults = ballKalman(p)
	x = zeros(1,p.nBug);
	realx = zeros(1,p.nBug);
	gam = p.gamma*ones(1,p.nBug);
	t = 0;

	innersteps = 1;
	outersteps = ceil(p.totalT);

	xmat = x;
	gammat = gam;
	realxmat = realx;
	tvec = t;

	k = 1;

	p.sigmaint = p.sigClean;
	while t < p.cleanTime
		[x gam realx t] = advanceBallKalman(p,x,gam,realx,t);
		k = k + 1;
		xmat(k,:) = x;
		gammat(k,:) = gam;
		realxmat(k,:) = realx;
		tvec(k) = t;
	end

	p.sigmaint = p.sigDirty;
	while t < p.cleanTime + p.dirtyTime
		[x gam realx t] = advanceBallKalman(p,x,gam,realx,t);
		k = k + 1;
		xmat(k,:) = x;
		gammat(k,:) = gam;
		realxmat(k,:) = realx;
		tvec(k) = t;
	end

	p.sigmaint = p.sigClean;
	while t < p.totalT
		[x gam realx t] = advanceBallKalman(p,x,gam,realx,t);
		k = k + 1;
		xmat(k,:) = x;
		gammat(k,:) = gam;
		realxmat(k,:) = realx;
		tvec(k) = t;
	end

	simresults.xmat = xmat;
	simresults.gammat = gammat;
	simresults.realxmat = realxmat;
	simresults.tvec = tvec;
end

function simresults = ballFixed(p)
	x = zeros(1,p.nBug);
	realx = zeros(1,p.nBug);
	t = 0;

	innersteps = ceil(0.05/p.dt);
	outersteps = ceil(p.totalT/(innersteps*p.dt));

	xmat = x;
	realxmat = realx;
	tvec = t;

	k = 1;

	p.sigmaint = p.sigClean;
	while t < p.cleanTime
		[x realx t] = advanceBallFixed(p,x,realx,t);
		k = k + 1;
		xmat(k,:) = x;
		realxmat(k,:) = realx;
		tvec(k) = t;
	end

	p.sigmaint = p.sigDirty;
	while t < p.cleanTime + p.dirtyTime
		[x realx t] = advanceBallFixed(p,x,realx,t);
		k = k + 1;
		xmat(k,:) = x;
		realxmat(k,:) = realx;
		tvec(k) = t;
	end

	p.sigmaint = p.sigClean;
	while t < p.totalT
		[x realx t] = advanceBallFixed(p,x,realx,t);
		k = k + 1;
		xmat(k,:) = x;
		realxmat(k,:) = realx;
		tvec(k) = t;
	end

	simresults.xmat = xmat;
	simresults.realxmat = realxmat;
	simresults.tvec = tvec;
end

%%%%%%%% Generic Functions
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

function stringy = howfar(pairs)
	stringy = '';
	for i = 1:length(pairs)
		stringy = [stringy num2str(pairs{i}(1)) ' of ' num2str(pairs{i}(2)) ' '];
	end
end

function makePretty()
	set(gca,'LineWidth',2)
	set(gca,'FontSize',20)
	% set(gca,'FontWeight','bold')
end