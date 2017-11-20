%HW7

% Problem 1: Modeling population growth
% The simplest model for a growing population assumes that each current
% individual has equal likelihood to divide, which yields a differential
% equation dx/dt = a*x where a is the division rate. It is easy to see that
% this will grow exponentially without bound. A simple modification is to
% assume that the growth rate slows done as the population reaches some
% maximum value N so that dx/dt = a*x*(1-x/N). Defining X = x/N, we have 
% dX/dt = a*X*(1-X).  
% Part 1. This equation has two fixed points at 0 and 1. Explain the
% meaning of these two points.

% The first fixed point, at X = 0, represents a population of 0, or 0% of
% the carrying capacity of the population. With 0 living individuals, no
% individuals are present to divide, so the population will have a growth
% rate of 0 and will remain stagnant at X = 0.
% The second fixed point, at X = 1, represents a population of N, or 100%
% of the carrying capacity of the population. With the population at the
% carrying capacity, the population cannot grow any more, so the growth
% rate of the population is 0, and the population will remain stagnant at
% X = 1.

% Part 2: Evaluate the stability of these fixed points. Does it depend on
% the value of the parameter a? 

% As explained above, each of these fixed points represents a steady state
% for the population. The fixed point at X = 0 is an unsteady state. If X
% rises above 0 at all, the rate dX/dt becomes positive, and the value of X
% continues to rise away from X = 0. The fixed point at X = 1 is a steady
% state. If X rises above 1, dX/dt becomes negative, and X decreases back
% towards 1. If X goes slightly below 1, dX/dt becomes positive, and X
% increases back towards 1.
% The value of the parameter a has no bearing on the stability of each of
% these fixed points, and does not factor into the above evaluation. a will
% only determine how quickly these small perturbations in X will correct
% towards the stable steady state.

% Part 3: Write a function that takes two inputs - the initial condition x0
% and the a parameter and integrates the equation forward in time. Make
% your code return two variables - the timecourse of X and the time
% required for the population to reach 99% of its maximum value. 

% See function written at bottom.

% Integrate and plot a sample result
[timecourse,t99] = pop_func(0.01, 5);
plot(timecourse(2,:), timecourse(1,:), 'r-')
title('Integrated Population Model, X_0 = 0.01, a = 5')
ylabel('Relative Population X = x/N')
xlabel('Time')
axis([0 5 0 1])
disp('For X0 = 0.01 and a = 5:')
disp(['t_99 = ' num2str(t99)])

% Part 4: Another possible model is to consider discrete generations
% instead allowing the population to vary continuously. e.g. X(t+1) = a*
% X(t)*(1-X(t)). Consider this model and vary the a parameter in the range 0
% < a <= 4. For each value of a choose 200 random starting points  0 < x0 < 1 
% and iterate the equation forward to steady state. For each final
% value Xf, plot the point in the plane (a,Xf) so that at the end you will
% have produced a bifucation diagram showing all possible final values of
% Xf at each value of a. Explain your results. 

% 0 = aX(t) - a(X(t)^2)
% ==>
% X(t) = (a + sqrt(a^2 - 0))/(2a) = 

figure
hold on
a = 0:4;
% Generate initial generations and initialize vector to hold final
% populations
x0s = rand(length(a), 200);
Xfs = zeros(size(x0s));

% Iterate through a values
for iter = 1:length(a)
    
    % Iterate through each of 200 sample start points
    for guess = 1:200
        Xf = x0s(iter, guess);
        
        % Take the population out to 100 generations to get steady state
        % value
        for gen = 1:100
            
            % Adjust population of new generation based on formula
            Xf = a(iter) * Xf * (1 - Xf);
            
        end
        
        % Store final population at each a value and starting condition
        Xfs(iter, guess) = Xf;
        
        % Plot this final value along with the iterated a
        plot(a(iter), Xf, 'r.', 'MarkerSize', 20)
        hold on
        
    end
    
end
xlabel('a')
ylabel('Xf')
title('Discrete Generational Population Growth, at Steady State')

% When analyzing the possible population values of the next generation, a
% steady state population value is found when:
% X(t) = aX(t) - aX(t)^2
% aX(t)^2 + (1-a)X(t) = 0
%
% When a = 0 or a = 1, steady state values for X(t) are only found at
% X(t) = 0.
% 0 * X(t)^2 + 1 * X(t) = 0  ==> X(t) = 0
% 1 * X(t)^2 + 0 * X(t) = 0  ==> X(t)^2 = 0
%
% When a = 2, X(t) reaches a steady state at X(t) = 0.5.
% 2 * X(t)^2 + (-1) * X(t) = 0  ==> 2X(t) = 1 ==> X(t) = 1/2
%
% When a = 3, X(t) reaches a steady state at X(t) = 0.667.
% 3 * X(t)^2 + (-2) * X(t) = 0  ==> 3X(t) = 2 ==> X(t) = 2/3
%
% When a = 4, X(t) reaches a steady state at X(t) = 0.75.
% 4 * X(t)^2 + (-3) * X(t) = 0  ==> 4X(t) = 3 ==> X(t) = 3/4
%
% However, with higher values of a, the difference between X(t) and X(t + dt)
% becomes greater. dX/dt = 2aX(t) + 1 - a takes on values further from 0
% for X(t) =/= X(t @ steady state), and X jumps around more. Once a = 4,
% the steady state is unstable, and population values randomly distribute
% themselves at many generations. X(t) will only stay at 0.75 if it lands
% on 0.75 exactly, which is highly unlikely.


% Problem 2. Genetic toggle switches. 
% Consider a genetic system of two genes A and B in which each gene
% product represses the expression of the other. Make the following
% assumptions: 
% a. Repression is cooperative:  each promotor region of one gene has 4
% binding sites for the other protein and that all of these need to be
% occupied for the gene to be repressed. 
% b. You can ignore the intermediate mRNA production so that the product of
% the synthesis of one gene can be assumed to directly repress the other
% c. the system is prefectly symmetric so that the degradation
% times, binding strengths etc are the same for both genes. 
% d. You can choose time and concentration scales so that all Michaelis
% binding constants and degradation times are equal to 1. 
%
% Part 1. Write down a two equation model (one for each gene product) for
% this system. Your model should have one free parameter corresponding to the
% maximum rate of expression of the gene, call it V. 
%
% dA/dt = (V + B^4)/(1 + B^4) - A
% dB/dt = (V + A^4)/(1 + A^4) - B
%
% Part 2. Write code to integrate your model in time and plot the results for V = 5 for two cases, 
% one in which A0 > B0 and one in which B0 > A0. 

AB0 = [2; 1];
timecourseA = gene_switch( AB0, 5 );
AB0 = [0.5; 4];
timecourseB = gene_switch( AB0, 5 );
figure
plot(timecourseA(3,:), timecourseA(1,:), 'r-')
hold on
plot(timecourseA(3,:), timecourseA(2,:), 'k-')
plot(timecourseB(3,:), timecourseB(1,:), 'r--')
plot(timecourseB(3,:), timecourseB(2,:), 'k--')
title('Integrated Model over Time')
ylabel('Gene Expression')
xlabel('Time')
legend('A, A0 > B0', 'B, A0 > B0', 'A, A0 < B0', 'B, A0 < B0')

% It is clear from the plotted integrated model that there are two steady
% state values for the Genes A and B to reach at V = 5. The factor that
% dictates which gene takes which steady state is the starting conditions.
% The gene that starts at a higher concentration will end up at a higher
% concentration at steady state, as it effectively repressed the other gene
% before the other gene can build up to the same expression.

% Part 3. By any means you want, write code to produce a bifurcation diagram showing all
% Iterate through a range of V values
for V = 0:0.05:8
    
    % Give a function for the steady state of gene B at a given V and B
    ss_b = @(B) V+((V+B^4)/(1+B^4))^4-B-B*((V+B^4)/(1+B^4))^4;
    
    % Iterate through a number of starting conditions for B to get all
    % roots present
    for Bo = 0:0.1:3
        
        % Use the fzero function to find the steady state root
        [rt,~,flag] = fzero(ss_b, [Bo]);
        if flag == 1
            % Plot the root on the bifurcation diagram
            plot(V,rt,'k.')
            hold on
        end
        
    end
    
end
title('Bifurcation Diagram')
ylabel('Fixed Points')
xlabel('V')

% The diagram shown here indicates that before about V = 3.5, there is only
% one steady state value for both Genes A and B. Their expressions are
% dictated more heavily by their repressions of each other, so neither is
% allowed to reach a larger expression without the other gene repressing it
% back down to its own level. 
% However, after about V = 3.5, basal expression becomes significant, and 
% two steady states arise, along with an intermediate unsteady state. The 
% unsteady state is when both expressions are exactly the same and
% repression and basal expression are evenly matched in each gene. The two
% steady states represent one gene dominating the other. If one gene has a
% higher expression than the other, the basal expression allows this one to
% reach a larger steady state and repress the other gene to a greater
% extent than that other gene represses it.
% This supports the data found in Part 2. At V = 5, two stable steady
% states were found.

% FUNCTIONS %

% Problem 1

function[ timecourse, t99 ] = pop_func( x0, a )

% Define integration parameters
dt = 0.01;
interval = [0 50];
nstep = (interval(2) - interval(1))/dt;
sol1(1) = x0;

% Iterate through all times
for ii = 2:nstep
    
    % Find new position for X
    sol1(ii) = sol1(ii-1) + dt * a * sol1(ii-1) * (1 - sol1(ii-1));
    
end

% Prepare data to return
times = linspace(interval(1), interval(2), nstep);
timecourse = [sol1; times];
t99 = times(find(timecourse > 0.99, 1));

end

% Problem 2

function[ timecourse ] = gene_switch( AB0, V )

% Define integration parameters
dt = 0.01;
interval = [0 20];
nstep = (interval(2) - interval(1))/dt;
sol1(:,1) = AB0;

% Iterate through all times
for ii = 2:nstep
    
    % Find new position for X
    sol1(:,ii) = sol1(:,ii-1) + dt * [((V + sol1(2,ii-1)^4)/(1 + sol1(2,ii-1)^4) - sol1(1,ii-1)); ((V + sol1(1,ii-1)^4)/(1 + sol1(1,ii-1)^4) - sol1(2,ii-1))];
    
end

% Prepare data to return
times = linspace(interval(1), interval(2), nstep);
timecourse = [sol1; times];

end
