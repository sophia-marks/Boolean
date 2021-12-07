% PatternFormation_StateLabels_AnyRule_MultipleScenarios.m
% Pattern formation plot for a network with N nodes governed by
% CA rules

clear;%clf;
N = 20; %#nodes
k = 3; %neighborhood size (all inputs including the node itself for ECA)
it = 20; %#iterations%
%Rules = [32 128 51 200 22 126 54 8 132 130 162 30 4 110 184 57 60 41 ]; 
Rules = 126;

marker = 2;
fsize = 14;
numlines = 1;
numcolumns = 1;
numplot = 1;
linewidth = 2;

%Here the parents and the initial condition are set the same for all rules
% generate the neighborhood for each node
parents_initial = parents_CA_includingthenode(N,k-1,N); %CA: k-1 adjacent nodes + the node itself (so a total of k)
%parents_initial = parents_random(N,k,N); %random choice of exactly k parents not including the node itself - so self-inputs are avoided here for obscure reasons

% generate the initial state of the network
A_initial = zeros(it,N); A_initial(1,floor(N/2))=1; %one black cell
%A_initial = ones(it,N);
%A_initial = randominitialstate(N); %random initial state

%This is where we focus on one rule at a time and create patterns and
%heatmaps
rulecount = 1;
for Rule = Rules %this will create twice as many figures so please be careful!
    clear rule_name parents A counts outputlabels
    rule_name = num2str(Rule);
    rule = binary_rule(rule_name, k-1);
    
    %reinitialize parents and A to be the same for all rules
    parents = parents_initial;
    A = A_initial;

% iterate the network and count the ON cells at each step
for j=2:it
[A(j,:), outputlabels(j,:)] = one_network_iteration_rules_StateLabels(A(j-1,:),parents,rule);
end

figure
A=A'; outputlabels = outputlabels';
[d1,d2] = size(A);
if d1 ~= N
    display('Trouble: d1 different from number of nodes');
end
if d2 ~= it
    display('Trouble: d2 different from number of iterations');
end
% plot the pattern formation and the density of ones
subplot(numlines,numcolumns,numplot);
for j=1:d1 %by node
    i = 1;
    if A(j,i) == 1
       plot(d1-(j-1),d2-(i-1),'sw','MarkerSize',marker,...
          'MarkerEdgeColor','k', 'MarkerSize',30);
      hold on;
    end
    for i=2:d2 %by iteration
        %pattern plot with labels
        if A(j,i) == 1
            plot(d1-(j-1),d2-(i-1),'sw','MarkerSize',marker,...
                'MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',30);
            hold on;
            a = cellstr(outputlabels(j,i));
            for m = 1:length(a{1,1})
                b(m) = str2num(a{1,1}(m));
            end
            plot(d1-(j-1),d2-(i-1),'sw','MarkerSize',marker,...
                 'MarkerEdgeColor',b,'MarkerFaceColor',b, 'MarkerSize',30);
            %plot(d1-(j-1),d2-(i-1),'sw','MarkerSize',marker,...
            %   'MarkerEdgeColor',str2num(outputlabels(j,i)),'MarkerFaceColor',str2num(outputlabels(j,i)), 'MarkerSize',30);
            text(d1-(j-1)-.3,d2-(i-1),outputlabels(j,i),'fontsize',8,'fontweight','bold');
            %text(d1-(j-1)-.2,d2-(i-1),outputlabels(j,i),'fontsize',8,'fontweight','bold', 'color',[rand,rand,rand]);
        end
    end %iteration loop
end %node loop
%hold off

%plot pattern
axis equal;
axis([1 it 1 N]);
axis off;
title(['Network evolution, N=',num2str(N), ', k=',num2str(k),...
    ', Rule ',rule_name], 'FontSize',14,'FontWeight','bold');

%prepare counts for the heatmap
for j = 1:d2 %counts by iterations
        counts(j,1)=sum(outputlabels(:,j)=="000");
        counts(j,2)=sum(outputlabels(:,j)=="001");
        counts(j,3)=sum(outputlabels(:,j)=="010");
        counts(j,4)=sum(outputlabels(:,j)=="011");
        counts(j,5)=sum(outputlabels(:,j)=="100");
        counts(j,6)=sum(outputlabels(:,j)=="101");
        counts(j,7)=sum(outputlabels(:,j)=="110");
        counts(j,8)=sum(outputlabels(:,j)=="111");
end
overallcount(rulecount,:) = sum(counts,1);
overallproportion(rulecount,:) = overallcount(rulecount,:)/sum(overallcount(rulecount,:));
rulecount = rulecount+1;

figure
%heatmap for proportion of states by iterations
%colormaps to be found at the bottom of the online page
%https://www.mathworks.com/help/matlab/ref/colormap.html#buc3wsn-1-map
%hsv turbo parula jet summer autumn .....
h = heatmap(counts/N, 'XData', ["000" "001" "010" "011" "100" "101" "110" "111"],'Colormap',jet);
end % rule for loop

figure
overallcount = [Rules' overallcount];
overallproportion = [Rules' overallproportion];
[dim1, dim2] = size(overallproportion);
for i = 1:dim1
    plot(overallproportion(i,2:dim2));
    hold on
end

%D = pdist(overallproportion(:,2:dim2)', 'euclidean');
   
T=clusterdata(overallcount,'Distance','euclidean',3);

M=[T';Rules];
figure(50)
scatter(Rules,T,'filled')
%display(M)

%Wolframclass=[1 1 2 2 3 3 4 1 1 2 2 3 4 4 2 2 3 4];

%Wuenscheclass=[1 1 1 1 1 1 1 2 1 2 2 2 1 2 3 3 3 3];

%RAND=rand_index(T',Wuenscheclass);

%display(RAND)

function ri = rand_index(p1, p2, varargin)
%RAND_INDEX Computes the rand index between two partitions.
%   RAND_INDEX(p1, p2) computes the rand index between partitions p1 and
%   p2. Both p1 and p2 must be specified as N-by-1 or 1-by-N vectors in
%   which each elements is an integer indicating which cluster the point
%   belongs to.
%
%   RAND_INDEX(p1, p2, 'adjusted') computes the adjusted rand index
%   between partitions p1 and p2. The adjustment accounts for chance
%   correlation.
    %% Parse the input and throw errors
    % Check inputs
    adj = 0;
    if nargin == 0
        error('Arguments must be supplied.');
    end
    if nargin == 1
        error('Two partitions must be supplied.');
    end
    if nargin > 3
        error('Too many input arguments');
    end
    if nargin == 3
        if strcmp(varargin{1}, 'adjusted')
            adj = 1;
        else
            error('%s is an unrecognized argument.', varargin{1});
        end
    end
    if length(p1)~=length(p2)
        error('Both partitions must contain the same number of points.');
    end
    % Check if arguments need to be flattened
    if length(p1)~=numel(p1)
        p1 = p1(:);
        warning('The first partition was flattened to a 1D vector.')
    end
    if length(p2)~=numel(p2)
        p2 = p2(:);
        warning('The second partition was flattened to a 1D vector.')
    end
    % Check for integers
    if isreal(p1) && all(rem(p1, 1)==0)
        % all is well
    else
        warning('The first partition contains non-integers, which may make the results meaningless. Attempting to continue calculations.');
    end
    if isreal(p2) && all(rem(p2, 1)==0)
        % all is well
    else
        warning('The second partition contains non-integers, which may make the results meaningless. Attempting to continue calculations.');
    end
	%% Preliminary computations and cleansing of the partitions
    N = length(p1);
    [~, ~, p1] = unique(p1);
    N1 = max(p1);
    [~, ~, p2] = unique(p2);
    N2 = max(p2);
    %% Create the matching matrix
    for i=1:1:N1
        for j=1:1:N2
            G1 = find(p1==i);
            G2 = find(p2==j);
            n(i,j) = length(intersect(G1,G2));
        end
    end
    %% If required, calculate the basic rand index
    if adj==0
        ss = sum(sum(n.^2));
        ss1 = sum(sum(n,1).^2);
        ss2 =sum(sum(n,2).^2);
        ri = (nchoosek2(N,2) + ss - 0.5*ss1 - 0.5*ss2)/nchoosek2(N,2);
    end
    %% Otherwise, calculate the adjusted rand index
    if adj==1
        ssm = 0;
        sm1 = 0;
        sm2 = 0;
        for i=1:1:N1
            for j=1:1:N2
                ssm = ssm + nchoosek2(n(i,j),2);
            end
        end
        temp = sum(n,2);
        for i=1:1:N1
            sm1 = sm1 + nchoosek2(temp(i),2);
        end
        temp = sum(n,1);
        for i=1:1:N2
            sm2 = sm2 + nchoosek2(temp(i),2);
        end
        NN = ssm - sm1*sm2/nchoosek2(N,2);
        DD = (sm1 + sm2)/2 - sm1*sm2/nchoosek2(N,2);
        % Special case to handle perfect partitions
        if (NN == 0 && DD==0)
            ri = 1;
        else
            ri = NN/DD;
        end
    end
    %% Special definition of n choose k
    function c = nchoosek2(a,b)
        if a>1
            c = nchoosek(a,b);
        else
            c = 0;
        end
    end
end



%scatter(overallproportion(:,1),overallproportion(:,2),100,T,'filled')
%title('Result of Clustering');
%scatter3(X(:,1),X(:,2),X(:,3),100,y,'filled')
%title('Randomly Generated Data in Three Clusters');


