% PROGRAM NAME: ps4huggett.m
clear, clc

% PARAMETERS
beta = .9932; %discount factor 
sigma = 1.5; % coefficient of risk aversion
b = 0.5; % replacement ratio (unemployment benefits)
y_s = [1, b]; % endowment in employment states
PI = [.97 .03; .5 .5]; % transition matrix


% ASSET VECTOR
a_lo = -2; %lower bound of grid points
a_hi = 5;%upper bound of grid points
num_a = 1000;

a = linspace(a_lo, a_hi, num_a); % asset (row) vector

% INITIAL GUESS FOR q
q_min = 0.98;
q_max = 1;
q_guess = (q_min + q_max) / 2;

% ITERATE OVER ASSET PRICES
aggsav = 1 ;
while abs(aggsav) >= 0.01 
    
    % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus, a', q_guess * a);
    cons = bsxfun(@plus, cons, permute(y_s, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret(cons<0)=-Inf;
    
    % INITIAL VALUE FUNCTION GUESS
    v_guessE = zeros(1, num_a);
    v_guessU = zeros(1, num_a);
    
    % VALUE FUNCTION ITERATION
    disE = 1; disU=1; tol = 1e-06;
    while  (disE > tol) && (disU > tol);
        
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
    expvalE = PI(1,1)*v_guessE + PI(1,2)*v_guessU;
    expvalU = PI(2,1)*v_guessE + PI(2,2)*v_guessU;
    
    value_matE = ret(:,:,1) + beta * repmat(expvalE, [num_a 1]);
    value_matU = ret(:,:,2) + beta * repmat(expvalU, [num_a 1]);
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        
    [vfn, pol_indxE] = max(value_matE, [], 2);
    vfnE = vfn';
    
    [vfn, pol_indxU] = max(value_matU, [], 2);
    vfnU = vfn';
  
    disE = max(abs(vfnE - v_guessE));
    disU = max(abs(vfnU - v_guessU));
    
    v_guessE = vfnE;
    v_guessU = vfnU;
    
    
    end;
    
    % KEEP DECSISION RULE
    pol_fnE = a(pol_indxE);
    pol_fnU = a(pol_indxU);
    % SET UP INITITAL DISTRIBUTION
        
    Mu=zeros(2,num_a); % 2 dimensions, one for each state
    Mu(:)=1/(2*num_a);
    
    % Combine Policy Indices to use with distribution Mu
    pol_indx = horzcat(pol_indxE, pol_indxU)';
    
    % ITERATE OVER DISTRIBUTIONS   
    distance=1; %Setup distance to check between previous and new distribution
    
    while distance >1e-06 
        
    [state_ind, a_ind, mass] = find(Mu > 0); % find non-zero indices
    
    MuNew = zeros(2, num_a);
    
    for j = 1:length(state_ind)
        
        apr_ind = pol_indx(state_ind(j), a_ind(j)); % which a prime does the policy fn prescribe?
       
        MuNew(:, apr_ind) = [PI(state_ind(j), :)'*Mu(state_ind(j),a_ind(j))];
        
        %MuNew(:, apr_ind) = MuNew(:, apr_ind) + ... % which mass of households goes to which exogenous state?
            %(PI(emp_ind(ii), :) * mass)';  
    end
    distance = max(max(abs(Mu-MuNew)));
    Mu=MuNew;
    end
    
    %Clear the Market
      %Market clear
   aggsav= Mu(1,:)*pol_fnE' + Mu(2,:)*pol_fnU';
   if aggsav>0;
        q_min=q_guess;
   else q_max=q_guess;
   end
   
end

figure(1)
 plot(a,vfnE,'blue')
 hold on
 plot(a,vfnU,'red')
 legend('Employed','Unemployed','location','southeast')
 title(['Value Function'])
 hold off
 
 figure(2)
 plot(a,pol_fnE,'blue')
 hold on
 plot(a,pol_fnU,'red')
 legend('Employed','Unemployed','location','southeast')
 title(['Policy Function'])
 hold off
 
% Wealth Distribution
wealth=[a+y_s(1);a+y_s(2)];
figure(3)
 bar(wealth(1,:),Mu(1,:),'blue')
 hold on
 bar(wealth(2,:),Mu(2,:),'red')
 xlim([-2 2.5]);
 legend('Employed','Unemployed','location','northwest')
 title('Distribution')
 hold off

 % Lorenz Curve and Gini
 
 wealth_n =reshape([a+y_s(1);a+y_s(2)],[2*num_a,1]);
 population = reshape(Mu',[2*num_a,1]);
 wealth_n(wealth_n<0)=0;
 
 figure(4)
 gini_w =gini(population, wealth_n,true);
 title(['Wealth Lorenz, Gini=',num2str(gini_w)])
 
 
 earnE= repmat(y_s(1),[1 num_a]);
 earnU= repmat(y_s(2),[1 num_a]);
 earnings= reshape([earnE ; earnU]',[2*num_a,1]);
 
 figure(5)
 gini_e=gini(population, earnings,true);
 title(['Earnings Lorenz, Gini=',num2str(gini_e)])
 
