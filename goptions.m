function OPTIONS=goptions(parain);
%GOPTIONS Default parameters used by the genetic algorithm GENETIC.
%
% Note that since the original version was written, the Matlab Optimization
% Toolbox now uses "optimset" to set generic optimization parameters, so
% this format is somewhat outdated.
% 
% The genetic algorithm parameters used for this implementation are:
%
%	OPTIONS(1)-Display flag:  0 = none, 1 = some, 2 = all  (Default: 1).
%	OPTIONS(2)-Termination bit string affinity value (Default: 0.90; set to zero to turn off)
%	OPTIONS(3)-Termination tolerance for fitness (Default: 0; not normally used).
%	OPTIONS(4)-Termination number of consecutive generations with same best
%	fitness (Default: 0; to use, set number, be sure OPTIONS(2) and OPTIONS(3) = 0).
%	OPTIONS(5)-Number of elite individuals (Default: 0; no elitism).
%	OPTIONS(6)-
%	OPTIONS(7)-
%	OPTIONS(8)-
%	OPTIONS(9)-
%	OPTIONS(10)-
% Genetic Algorithm-specific inputs
%	OPTIONS(11)-Population size (fixed)
%	OPTIONS(12)-Probability of crossover
%	OPTIONS(13)-Probability of mutation
%	OPTIONS(14)-Maximum number of generations, always used as safeguard
%	(Default: 200).
%	
%
% Explanation of defaults:
%	The default algorithm displays statistical information for each
%	generation by setting OPTIONS(1) = 1.  Plots are produced when
%	OPTIONS(1) = 2.  
%   The OPTIONS(2) flag is originally set for termination criterion based
%   on X; here it is used if bit string affinity is selected.
%   The default fitness function termination tolerance,
%	OPTIONS(3), is set to 0, which terminates the optimization when 5
%	consecutive best generation fitness values are the same.  A positive
%	value terminates the optimization when the normalized difference
%	between the previous fitness and current generation fitness is less
%	than the tolerance.  See the code for details.
%  OPTIONS(4) has a default value of 5; this means if the best fitness
%	value in the population is unchanged for 5 consecutive generations
%  the GA is terminated.
%       The default algorithm uses a fixed population size, OPTIONS(11),
%       and no generational overlap.  The default population size is 30.
%	Three genetic operations:  selection, crossover, and mutation are
%	used for procreation.
%	The default selection scheme is tournament selection.
%       Crossover occurs with probability Pc=OPTIONS(12).  The default
%	crossover scheme is uniform crossover with Pc = 0.5.
%	Each allele of the offspring mutates independently with probability
%       Pm=OPTIONS(13); here the default is 0.01.
%       The default number of maximum generations, OPTIONS(14) is 200.
%
%   Last modified by Bill Crossley 07/20/04

% The following lines have been commented out by Steven Lamberson.
% They have been changed to what is seen below them. (06/30/06).
% This change was made in order to fix the following problems:
%   1 - code changed user supplied options(1)=0 to options(1)=1
%   2 - code changed user supplied options(2)=0 to options(2)=0.9

%if nargin<1; parain = []; end
%sizep=length(parain);
%OPTIONS=zeros(1,14);
%OPTIONS(1:sizep)=parain(1:sizep);
%default_options=[1,0.9,0,0,0,0,0,0,0,0,0,0,0,200];
%OPTIONS=OPTIONS+(OPTIONS==0).*default_options

if nargin<1; parain = []; end
sizep=length(parain);
OPTIONS=zeros(1,14)-1;
OPTIONS(1:sizep)=parain(1:sizep);
default_options=[1,0.9,0,0,0,0,0,0,0,0,0,0,0,200];
for i = 1:length(OPTIONS)
    if OPTIONS(i) == -1
        OPTIONS(i) = default_options(i);
    end
end

