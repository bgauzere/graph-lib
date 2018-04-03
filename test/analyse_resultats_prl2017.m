datasets={'acyclic','mao'}
costs={'metric','nonmetric'}
solvers={'FBP','FBP0','SFBP','EBP','ECBP','FLWC'}
M=length(solvers)
results_lb=zeros(2,2,M);
results_ub=zeros(2,2,M);
results_lb_times=zeros(2,2,M);
results_ub_times=zeros(2,2,M);
%results pour la lower bound
for d = 1:length(datasets)
    for c =1:length( costs)
        for solver = 1:length(solvers)
            prefix=[datasets{d} '_' costs{c} '_' solvers{solver}];
            lb=load([prefix '_lb.mat'],'-ascii');
            results_lb(d,c,solver) = mean(mean(lb));
            
            times=load([prefix '_lb_times.mat'],'-ascii');
            results_lb_times(d,c,solver) = mean(mean(times));
        end
    end
end

%upper bounds
for d = 1:length(datasets)
    for c =1:length( costs)
        for solver = 1:length(solvers)
            prefix=[datasets{d} '_' costs{c} '_' solvers{solver}];
            ub=load([prefix '_ub.mat'],'-ascii');
            results_ub(d,c,solver) = mean(mean(ub));
            
            times=load([prefix '_ub_times.mat'],'-ascii');
            results_ub_times(d,c,solver) = mean(mean(times));
            
        end
    end
end

%defaults:
defaults_lb=zeros(2,M);
defaults_ub=zeros(2,M);
c=2
for d = 1:length(datasets)
    exact = load([datasets{d} '_exact'],'-ascii');
    N=sqrt(size(exact,1));
    exact=reshape(exact,N,N);    
    for solver = 1:length(solvers)
        prefix=[datasets{d} '_' costs{c} '_' solvers{solver}];
        lb=load([prefix '_lb.mat'],'-ascii');
        defaults_lb(d,c,solver) = sum(sum(lb > exact));
        
        
        prefix=[datasets{d} '_' costs{c} '_' solvers{solver}];
        ub=load([prefix '_ub.mat'],'-ascii');
        defaults_ub(d,c,solver) = sum(sum(ub < exact));
    end
end

fprintf('defaults LB Acyclic:\n')
for solver = 1:length(solvers)
    fprintf('%s : %d \n',solvers{solver}, defaults_lb(1,2,solver));
end


fprintf('defaults UB Acyclic:\n')
for solver = 1:length(solvers)
    fprintf('%s : %d \n',solvers{solver}, defaults_ub(1,2,solver));
end


fprintf('defaults LB MAO:\n')

for solver = 1:length(solvers)
    fprintf('%s : %d \n',solvers{solver}, defaults_lb(2,2,solver));
end


fprintf('defaults UB MAO:\n')

for solver = 1:length(solvers)
    fprintf('%s : %d \n',solvers{solver}, defaults_ub(2,2,solver));
end

fprintf('Export Latex Table 4\n')
for solver = 1:length(solvers)
    fprintf('%s & %e & %d & %e & %d \n',solvers{solver}, ...
            results_lb_times(1,2,solver), defaults_lb(1,2,solver), ...
            results_lb_times(2,2,solver), defaults_lb(2,2,solver))
end



fprintf('Export Table UB \n')
fprintf('Metric cost \n')

for solver = 1:length(solvers)
    fprintf('%s & %.2e & %.2f & %.2e & %.2f \n',solvers{solver}, ...
            results_ub_times(1,1,solver), results_ub(1,1,solver), ...
            results_ub_times(2,1,solver), results_ub(2,1,solver))
end

fprintf('Non Metric cost \n')
for solver = 1:length(solvers)
    
    fprintf('%s & %.2e & %.2f & %.2e & %.2f \n',solvers{solver}, ...
            results_ub_times(1,2,solver), results_ub(1,2,solver), ...
            results_ub_times(2,2,solver), results_ub(2,2,solver))
end
