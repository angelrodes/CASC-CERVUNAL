function [mu,sigma] = BFG(x,y)
% BFG calculates the Best Gaussian Fit of a probability distribution
% defined by x and y, where y are the probabilities.
%
% Angel Rodes, 2020
% www.angelrodes.com

% About the function name: The acronym for "Best Gaussian Fit" should be "BGF". However, when I created the function file, I made an unconscious mistake, probably influenced by my Best Friends in Glasgow, or the Big Friendly Giant that leads the SUERC-cosmo group. For that reason (and also because I referenced this function in many other scripts), I decided to leave the "mistake" as it is :)

% conver to array
x=x(:); y=y(:);

% check length
if size(x)~=size(y)
    warning('BGF: x and y should be the same length.')
end

% define gussian pdf
normalprob=@( x, mu, s )...
    1/(s*(2*pi)^0.5)*...
    exp(-(x-mu).^2/(2*s^2));

% define R^2
rsquare=@(mu,s)1-(sum(abs(y-normalprob( x, mu, s)*sum(y)/sum(normalprob( x, mu, s))).^2))/sum((mean(y)-y).^2);

% define models
nmodels=10000; % max number of models to run
mi=zeros(1,nmodels);
si=zeros(1,nmodels);
ri=zeros(1,nmodels)-9999;

% define initial limits
MUs=[min(x),max(x)];
SDs=[min(diff(x)),max(x)-min(x)];

% define tolerance
tol=0.001;

% run models
h = waitbar(0,'BGF');
n=0;
convergence=0;
while convergence==0
    n=n+1;
    if mod(n,100)==0
        waitbar(n/nmodels,h)
    end
    mi(n)=rand*(MUs(2)-MUs(1))+MUs(1);
    si(n)=exp(rand*(log(SDs(2))-log(SDs(1)))+log(SDs(1)));
    ri(n)=rsquare(mi(n),si(n));
    if mod(n,100)==0 && n>1000
        sortedri=sort(ri(1:n),'descend');
        limitvalue=sortedri(round(n/10)); % converge to best 10% models
        selec=(ri>=limitvalue);
        MUs=[min(mi(selec)),max(mi(selec))];
        SDs=[min(si(selec)),max(si(selec))];
            if std(MUs)/mean(MUs)<tol && std(SDs)/mean(SDs)<tol
                convergence=1;
            end
    end
    if n==nmodels
        convergence=1;
    end
end
close(h)
% selec=find(ri==max(ri));
% selec=selec(end);
% mu=mi(selec);
% sigma=si(selec);
limitvalue=sortedri(100); % best 100 models (avoid random artifacts due to x sampling limitations)
selec=(ri>=limitvalue);
mu=mean(mi(selec));
sigma=mean(si(selec));
end

