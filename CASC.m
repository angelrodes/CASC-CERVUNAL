%% Cosmo-Ages Sequence Calculator
%  A simple bayesian approach for a set of gaussian ages in a sequence.
%
% Age probabilities for each unit (e.g. each morainte) are calculated as:
% P(t|T+>T>T-) ~ P(T|t) * P(t|t<T+) * P(t|t>T-)
% Where T+, T and T- are the measured data corresponding to the older,
% contemporary, and younger units.
% P(t|t<T+) and P(t|t>T-) are the product of
% old-to-young and young-to-old cumulative sum of older and younger 
% age distributions from other units, scaled between 0 and 1.
% All distributions are scaled to fit a total probability of 1 for each
% unit.
% Gaussians are fitted to resulting distributions to produce symmetric ages
% for the units.
%
% Angel Rodes
% SUERC 2020



%% clear previous data and plots
clear
close all hidden
clc

%% Version
script_version='0.4.1';

%% Import data
fid = fopen('CERVUNAL_internal.csv');
% Sample name , age , uncertainty , moraine name ,  seqential order
mydata = textscan(fid, '%s %f %f %s %f',...
    'HeaderLines', 1,'Delimiter',',');
fclose(fid);

sample_names=mydata{1};
ages=mydata{2};
uncertainties=mydata{3};
unit_names=mydata{4};
unit_order=mydata{5};

%% Minimum time between units
% % Ask minimum time between sequential units
% % Ask for offset. Default value = 0
% answer=inputdlg('Minimum time between units:',...
%     'Same units as in data.csv!',[1 50],{'0'});
% % Convert answer to number
% minimum_time_between_units=str2num(answer{1});
% % Check that this is a number
% if isempty(minimum_time_between_units)
%     error('Minimum time between units is not a number!')
% end

minimum_time_between_units=0;

%% Define useful functions

% Normal function
normal_probs=@(mu,sigma,x)...
    1/(sigma*(2*pi())^0.5)*exp(-(x-mu).^2/(2*sigma^2));

% One-sigma-values function
one_sigma_percent=0.682689492;
one_sigma_threshold=@(P)...
    interp1(sort(P,'descend'),...
    find(cumsum(sort(P,'descend')/sum(P))>=one_sigma_percent,1,'first'));
one_sigma_selection=@(P)...
    P>=one_sigma_threshold(P);


%% Define t space
% Maximum number of units
max_units=max(unit_order)-min(unit_order)+2;
% Minimum and maximum age values
min_t=min(ages-2*uncertainties-minimum_time_between_units*max_units);
max_t=max(ages+2*uncertainties+minimum_time_between_units*max_units);
% Time resolution
minimum_n_data=1000;
maximum_n_data=1e6;
minimum_points_per_age=10;
t_step=max(...
    min(min(2*uncertainties)/minimum_points_per_age,...
    (max_t-min_t)/minimum_n_data),...
    (max_t-min_t)/minimum_n_data);
% Time space
t=min_t:t_step:max_t;

%% Define initial probability matrix
unique_units=[min(unit_order)-1:max(unit_order)+1]';
% In principle, all ages are equally probable
%  and the sum of probabilities for each unit is 1
KDE=ones(length(unique_units),length(t))/length(t);
% KDE is for kernel density estimation

% for each unit with data
for seq=unique(unit_order)'
    % select samples
    sel=find(seq==unit_order);
    if ~isempty(sel)
        % calulate the kernel density estimation of the unit
        kdei=t.*0;
        for sample=sel'
            kdei=kdei+normal_probs(ages(sample),uncertainties(sample),t);
        end
        KDE((seq==unique_units),:)=kdei/sum(kdei);
        % Scaling to sum(kdei) forces the sum of each distribution = 1
    end
end

%% Calculate sequence kernel density estimations (SKDE)
% start a new matrix
SKDE=ones(size(KDE))/length(t);

% for each unit
for seq=unique_units'
    % get the probabilities from older units
    P_older=1+0.*t;
    for seq2=unique_units(unique_units<seq)'
        kdei=KDE((seq2==unique_units),:);
        offset=(seq-seq2)*minimum_time_between_units;
        kdei_with_offset=interp1(t-offset,kdei,t,'nearest','extrap');
        kdeolder=kdei_with_offset;
        P_older=P_older.*(1-cumsum(kdeolder/sum(kdeolder)));
    end
    
    
    % get the probabilities from younger units
    P_younger=1+0.*t;
    for seq2=unique_units(unique_units>seq)'
        kdei=KDE((seq2==unique_units),:);
        offset=(seq2-seq)*minimum_time_between_units;
        kdei_with_offset=interp1(t+offset,kdei,t,'nearest','extrap');
        kdeyounger=kdei_with_offset;
        P_younger=P_younger.*cumsum(kdeyounger/sum(kdeyounger));
    end
    
    kdei=KDE((seq==unique_units),:);
    
    skdei=kdei.*P_younger.*P_older;

    
    if sum(skdei)>0 % ignore unit if no-solution
        SKDE((seq==unique_units),:)=skdei/sum(skdei);
        % Scaling to sum(skdei) forces the sum of each distribution = 1
    end
end

%% Calculate ranges and display results
% Set precision
significant_decimal=-floor(log10(min(uncertainties)))+1;
precision_multiplier=10^(significant_decimal);

% Define min and max t limits for plotting
min_tplot=min(ages-1*uncertainties);
max_tplot=max(ages+1*uncertainties);


% Credits
disp(['Cosmo-Ages Sequence Calculator v.' script_version])
disp('Angel Rodes, SUERC 2020')
disp('--')

% Results
if minimum_time_between_units~=0
    disp(['Minimum time between units: '...
        num2str(minimum_time_between_units)])
end
disp('One-sigma ranges:')
for seq=sort(unique_units)'
    if seq==min(unique_units)
        disp(['  Unit ' num2str(seq) ': Start'])
    elseif seq==max(unique_units)
        disp(['  Unit ' num2str(seq) ': End'])
    elseif sum(seq==unit_order)>0
        name_seq=unit_names{find(seq==unit_order,1,'first')};
        disp(['  Unit ' num2str(seq) ': ' name_seq])
    else
        disp(['  Unit ' num2str(seq)])
    end
    
    kdei=KDE((seq==unique_units),:);
    sel=one_sigma_selection(kdei);
    if min(t(sel))==min(t)
        minrange=-Inf;
    else
        minrange=min(t(sel));
    end
    if max(t(sel))==max(t)
        maxrange=Inf;
    else
        maxrange=max(t(sel));
    end
    disp(['    Input data: [ '...
        num2str(round(minrange*precision_multiplier)/precision_multiplier)....
        ' - '...
        num2str(round(maxrange*precision_multiplier)/precision_multiplier)...
        ' ]'])

    kdei=SKDE((seq==unique_units),:);
    sel=one_sigma_selection(kdei);
    if min(t(sel))==min(t)
        minrange=-Inf;
    else
        minrange=min(t(sel));
    end
    if max(t(sel))==max(t)
        maxrange=Inf;
    else
        maxrange=max(t(sel));
    end
    disp(['      Sequence: [ '...
        num2str(round(minrange*precision_multiplier)/precision_multiplier)....
        ' - '...
        num2str(round(maxrange*precision_multiplier)/precision_multiplier)...
        ' ]'])
    % Calculte BGF
    [mu,sigma] = BFG(t,kdei);
    Bgf(find(unique_units==seq))=mu;
    dBgf(find(unique_units==seq))=sigma;
    significant_decimal_BGF=-floor(log10(sigma))+1;
    precision_multiplier_BGF=10^(significant_decimal_BGF);
    disp(['      Sequence BGF: '...
        num2str(round(mu*precision_multiplier_BGF)/precision_multiplier_BGF)....
        ' ± '...
        num2str(round(sigma*precision_multiplier_BGF)/precision_multiplier_BGF)...
        ' '])
    
    % Include Start and End limits in the plots
    min_tplot=min(min_tplot,maxrange-median(uncertainties));
    max_tplot=max(max_tplot,minrange+median(uncertainties));
    
end

%% Plot camelplots
% Scale probabiltiies so max(P)=0.8
% for plotting purposes
scale_kdes=0.8/max(max(KDE(:)),max(SKDE(:)));

% define minimum and maximum plot limits
% min_tplot=min(ages-3*uncertainties);
% max_tplot=max(ages+3*uncertainties);
min_seqplot=min(unique_units)+1;
max_seqplot=max(unique_units)-1;

figure
hold on

% for each unit
for seq=unique_units'
    
    selsamples=(seq==unit_order);
    seqpos=find(unique_units==seq);

   

    % plot initial KDEs
    x=t;
    kdei=KDE((seq==unique_units),:);
    if length(unique(kdei))>1 % plot only if significant data
        y=seq+kdei*scale_kdes;
        plot(x,y,'-','LineWidth',2,'Color',[0.4 0.4 0.4])
        % plot one-sigma bars
        sel=one_sigma_selection(kdei);
        vertical_offset=0.05;
        %         plot(x.*1./sel,(seq-vertical_offset)+y.*0*1./sel,'-b')
        yline=(seq-vertical_offset)+y.*0*1./sel;
        yline=mean(yline(~isnan(yline)));
%         plot([Bgf(seqpos)-dBgf(seqpos),Bgf(seqpos)+dBgf(seqpos)],yline*[1,1],'-k')
%         text(Bgf(seqpos)-dBgf(seqpos),yline,'<','HorizontalAlignment','center','VerticalAlignment','middle')
%         text(Bgf(seqpos)+dBgf(seqpos),yline,'>','HorizontalAlignment','center','VerticalAlignment','middle')
        quiver(Bgf(seqpos),yline,dBgf(seqpos),0,'-k')
        quiver(Bgf(seqpos),yline,-dBgf(seqpos),0,'-k')
        name_seq=unit_names(find(seq==unit_order,1,'first'));
        text(Bgf(seqpos),yline,name_seq,'Color','k',...
            'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','w')
    end

     if sum(selsamples)>0 % plot only if significant data
        % plot samples
        n=0;
        for sample=find(selsamples)'
            yi=normal_probs(ages(sample),uncertainties(sample),t);
            scale_this=1/sum(selsamples);
            yiplot=seq+yi*scale_this;
            plot(t,yiplot,'-k')
            text(ages(sample), max(yiplot),[' ' sample_names(sample)],...
                'HorizontalAlignment','left', 'VerticalAlignment', 'middle', 'Rotation', 90,...
                'FontSize',5)
%             plot([ages(sample)-uncertainties(sample),...
%                 ages(sample)+uncertainties(sample)],...
%                 [y,y],'-k')
%             plot(ages(sample),y,'xk')
        end
    end
    
    
    % plot bottom line
    plot(t,seq+t.*0,'-','Color',[0.8 0.8 0.8],'LineWidth',2)
end


xlim([min_tplot max_tplot])
ylim([min_seqplot-0.2 max_seqplot+0.85])
% ylabel('Sequence')
xlabel('Age (ka)')
set(gca,'YTick',[])
if minimum_time_between_units~=0
    title(['Minimum time between units: '...
        num2str(minimum_time_between_units)])
end
box on
grid on

