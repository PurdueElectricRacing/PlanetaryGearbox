function [xopt,fopt,stats,nfit,fgen,lgen,lfit] = GA550(fun, ...
    x0,options,vlb,vub,bits,P1,P2,P3,P4,P5,P6,P7P,P8,P9,P10)
%GA550 minimizes a fitness function using a simple genetic algorithm.
%
%	X=GA550('FUN',X0,OPTIONS,VLB,VUB) uses a simple  
%       genetic algorithm to find a minimum of the fitness function 
%       FUN.  FUN can be a user-defined M-file: FUN.M, or it can be a 
%	string containing the function itself.  The user may define all
%       or part of an initial population X0. Any undefined individuals 
%	will be randomly generated between the lower and upper bounds
%	(VLB and VUB).  If X0 is an empty matrix, the entire initial
%	population will be randomly generated.  Use OPTIONS to specify 
%	flags, tolerances, and input parameters.  Type HELP GOPTIONS
%       for more information and default values.
%
%	X=GA550('FUN',X0,OPTIONS,VLB,VUB,BITS) allows the user to 
%	define the number of BITS used to code non-binary parameters
%	as binary strings.  Note: length(BITS) must equal length(VLB)
%	and length(VUB).  If BITS is not specified, as in the previous 
%	call, the algorithm assumes that the fitness function is 
%	operating on a binary population.
%
%	X=GA550('FUN',X0,OPTIONS,VLB,VUB,BITS,P1,P2,...) allows up 
%	to ten arguments, P1,P2,... to be passed directly to FUN.
%	F=FUN(X,P1,P2,...). If P1,P2,... are not defined, F=FUN(X).
%
%	[X,FOPT,STATS,NFIT,FGEN,LGEN,LFIT]=GA550(<ARGS>)
%          X       - design variables of best ever individual
%          FOPT    - fitness value of best ever individual
%          STATS   - [min mean max stopping_criterion] fitness values 
%                    for each generation
%          NFIT	 - number of fitness function evalations
%          FGEN    - first generation population
%          LGEN    - last generation population
%          LFIT    - last generation fitness
%
%       The algorithm implemented here is based on the book: Genetic
%       Algorithms in Search, Optimization, and Machine Learning,
%       David E. Goldberg, Addison-Wiley Publishing Company, Inc.,
%       1989.
%
%	Originally created on 1/10/93 by Andrew Potvin, Mathworks, Inc. 
%	Modified on 2/3/96 by Joel Grasmeyer.
%   Modified on 11/12/02 by Bill Crossley.
%   Modified on 7/20/04 by Bill Crossley.

% Make best_feas global for stopping criteria (4/13/96)
global best_feas
global gen
global fit_hist
% Load input arguments and check for errors
if nargin<4,
    error('No population bounds given.')
elseif (size(vlb,1)~=1) | (size(vub,1)~=1),
    % Remark: this will change if algorithm accomodates matrix variables
    error('VLB and VUB must be row vectors')
elseif (size(vlb,2)~=size(vub,2)),
    error('VLB and VUB must have the same number of columns.')
elseif (size(vub,2)~=size(x0,2)) & (size(x0,1)>0),
    error('X0 must all have the same number of columns as VLB and VUB.')
elseif any(vlb>vub),
    error('Some lower bounds greater than upper bounds')
else
    x0_row = size(x0,1);
    for i=1:x0_row,
        if any(x0(x0_row,:)<vlb) | any(x0(x0_row,:)>vub),
            error('Some initial population not within bounds.')
        end % if initial pop not within bounds
    end % for initial pop
end % if nargin<4   

if nargin<6,
    bits = [];
elseif (size(bits,1)~=1) | (size(bits,2)~=size(vlb,2)),
    % Remark: this will change if algorithm accomodates matrix variables
    error('BITS must have one row and length(VLB) columns')
elseif any(bits~=round(bits)) | any(bits<1),
    error('BITS must be a vector of integers >0')
end % if nargin<6

% Form string to call for function evaluation
if ~( any(fun<48) | any(fun>122) | any((fun>90) & (fun<97)) | ...
        any((fun>57) & (fun<65)) ), 
    % Only alphanumeric characters implies that 'fun' is a separate m-file
    evalstr = [fun '(x'];
    for i=1:nargin-6,
        evalstr = [evalstr,',P',int2str(i)];
    end
else
    % Non-alphanumeric characters implies that the function is contained 
    % within the single quotes
    evalstr = ['(',fun];
end

% Determine all options
% Remark: add another options index for type of termination criterion
if size(options,1)>1,
    error('OPTIONS must be a row vector')
else
    % Use default options for those that were not passed in
    options = goptions(options);
end
PRINTING = options(1);
BSA = options(2);
fit_tol = options(3);
nsame = options(4)-1;
elite = options(5);

% Since operators are tournament selection and uniform crossover and
% default coding is Gray / binary, set crossover rate to 0.50 and use
% population size and mutation rate based on Williams, E. A., and Crossley,
% W. A., "Empirically-derived population size and mutation rate guidelines
% for a genetic algorithm with uniform crossover," Soft Computing in
% Engineering Design and Manufacturing, 1998.  If user has entered values
% for these options, then user input values are used.
if options(11) == 0,
    pop_size = sum(bits) * 4;
else
    pop_size = options(11);
end
if options(12) == 0,
    Pc = 0.5;
else
    Pc = options(12);
end
if options(13) == 0,
    Pm = (sum(bits) + 1) / (2 * pop_size * sum(bits));
else
    Pm = options(13);
end
max_gen = options(14);
% Ensure valid options: e.q. Pc,Pm,pop_size,max_gen>0, Pc,Pm<1
if any([Pc Pm pop_size max_gen]<0) | any([Pc Pm]>1),
    error('Some Pc,Pm,pop_size,max_gen<0 or Pc,Pm>1')
end

% Encode fitness (cost) function if necessary
ENCODED = any(any(([vlb; vub; x0]~=0) & ([vlb; vub; x0]~=1))) |  ....
    ~isempty(bits);
if ENCODED,
    [fgen,lchrom] = encode(x0,vlb,vub,bits);
else
    fgen = x0;
    lchrom = size(vlb,2);
end

% Display warning if initial population size is odd
if rem(pop_size,2)==1,
    disp('Warning: Population size should be even.  Adding 1 to population.')
    pop_size = pop_size +1;
end

% Form random initial population if not enough supplied by user
if size(fgen,1)<pop_size,
    fgen = [fgen; (rand(pop_size-size(fgen,1),lchrom)<0.5)];
end
xopt = vlb;
nfit = 0;
new_gen = fgen;
isame = 0;
bitlocavg = mean(fgen,1);  % initial bit string affinity
BSA_pop = 2 * mean(abs(bitlocavg - 0.5));
fopt = Inf;
stats = [];

% Header display
if PRINTING>=1,
    if ENCODED,
        disp('Variable coding as binary chromosomes successful.')
        disp('')
        fgen = decode(fgen,vlb,vub,bits);
    end
    disp('                   Fitness statistics')
    if nsame > 0
        disp('Generation Minimum      Mean         Maximum       isame')
    elseif BSA > 0
        disp('Generation Minimum      Mean         Maximum       BSA')
    else
        disp('Generation Minimum      Mean         Maximum       not used')
    end
end

% Set up main loop
STOP_FLAG = 0;
for generation = 1:max_gen+1,
    old_gen = new_gen;
    
    % Decode binary strings if necessary
    if ENCODED,
        x_pop = decode(old_gen,vlb,vub,bits);
    else
        x_pop = old_gen;
    end
    
    % Get fitness of each string in population
    for i = 1:pop_size,
        x = x_pop(i,:);
        fitness(i) = eval([evalstr,')']);
        nfit = nfit + 1;
    end
    
    % Store minimum fitness value from previous generation (except for
    % initial generation)
    if generation > 1,
        min_fit_prev = min_fit;
        min_gen_prev = min_gen;
        min_x_prev = min_x;
    end
    
    % identify worst (maximum) fitness individual in current generation
    [max_fit,max_index] = max(fitness);
    
    % impose elitism - currently only one individual; this replaces worst
    % individual of current generation with best of previous generation
    if (generation > 1 & elite > 0),   
        old_gen(max_index,:) = min_gen_prev;
        x_pop(max_index,:) = min_x_prev;
        fitness(max_index) = min_fit_prev;
    end
     
    % identify best (minimum) fitness individual in current generation and
    % store bit string and x values
    [min_fit,min_index] = min(fitness);
    min_gen = old_gen(min_index,:);
    min_x = x_pop(min_index,:);
    
    % Store best fitness and x values
    if min_fit < fopt,
        fopt = min_fit;
        xopt = min_x;
    end
    
    % Compute values for isame or BSA_pop stopping criteria
    if nsame > 0
        if generation > 1
            if min_fit_prev == min_fit
                isame = isame + 1;
            else
                isame = 0;
            end
        end
    elseif BSA > 0
        bitlocavg = mean(old_gen,1);
        BSA_pop = 2 * mean(abs(bitlocavg - 0.5));
    end
    
    
    % Calculate generation statistics
    if nsame > 0
        stats = [stats; generation-1,min(fitness),mean(fitness), ...
            max(fitness), isame];
    elseif BSA > 0
        stats = [stats; generation-1,min(fitness),mean(fitness), ...
            max(fitness), BSA_pop];
    else
        stats = [stats; generation-1,min(fitness),mean(fitness), ...
            max(fitness), 0];
    end
    
    % Display if necessary
    if PRINTING>=1,
        disp([sprintf('%5.0f %12.6g %12.6g %12.6g %12.6g', stats(generation,1), ...
                stats(generation,2),stats(generation,3), stats(generation,4),...
                stats(generation,5))]);
    end
    
    % Check for termination
    % The default termination criterion is bit string affinity.  Also
    % available are fitness tolerance across five generations and number of
    % consecutive generations with same best fitness.  These can be used
    % concurrently.
    if fit_tol>0,    % if fit_tol > 0, then fitness tolerance criterion used
        if generation>5,
            % Check for normalized difference in fitness minimums
            if stats(generation,1) ~= 0,
                if abs(stats(generation-5,1)-stats(generation,1))/ ...
                        stats(generation,1) < fit_tol
                    if PRINTING >= 1
                        fprintf('\n')
                        disp('GA converged based on difference in fitness minimums.')
                    end
                    lfit = fitness;
                    if ENCODED,
                        lgen = x_pop;
                    else
                        lgen = old_gen;
                    end
                    return
                end
            else
                if abs(stats(generation-5,1)-stats(generation,1)) < fit_tol
                    if PRINTING >= 1
                        fprintf('\n')
                        disp('GA converged based on difference in fitness minimums.')
                    end
                    lfit = fitness;
                    if ENCODED,
                        lgen = x_pop;
                    else
                        lgen = old_gen;
                    end
                    return
                end
            end
        end
    elseif nsame > 0,    % consecutive minimum fitness value criterion
            if isame == nsame
                if PRINTING >= 1
                    fprintf('\n')
                    disp('GA stopped based on consecutive minimum fitness values.')
                end
                lfit = fitness;
                if ENCODED,
                    lgen = x_pop;
                else
                    lgen = old_gen;
                end
                return
            end
    elseif BSA > 0,  % bit string affinity criterion
        if BSA_pop >= BSA,
            if PRINTING >=1
                fprintf('\n')
                disp('GA stopped based on bit string affinity value.')
            end
            lfit = fitness;
            if ENCODED,
                lgen = x_pop;
            else
                lgen = old_gen;
            end
            return
        end
    end
    
    % Tournament selection
    new_gen = tourney(old_gen,fitness);
    
    % Crossover
    new_gen = uniformx(new_gen,Pc);
    
    % Mutation
    new_gen = mutate(new_gen,Pm);
    
    % Always save last generation.  This allows user to cancel and
    % restart with x0 = lgen
    if ENCODED,
        lgen = x_pop;
    else
        lgen = old_gen;
    end
    
    
end % for max_gen

% Maximum number of generations reached without termination
lfit = fitness;
if PRINTING>=1,
    fprintf('\n')
    disp('Maximum number of generations reached without termination')
    disp('criterion met.  Either increase maximum generations')
    disp('or ease termination criterion.')
end


% end genetic

function [gen,lchrom,coarse,nround] = encode(x,vlb,vub,bits)
%ENCODE Converts from variable to binary representation.
%	[GEN,LCHROM,COARSE,nround] = ENCODE(X,VLB,VUB,BITS)
%       encodes non-binary variables of X to binary.  The variables
%       in the i'th column of X will be encoded by BITS(i) bits.  VLB
%       and VUB are the lower and upper bounds on X.  GEN is the binary
%       representation of these X.  LCHROM=SUM(BITS) is the length of
%       the binary chromosome.  COARSE(i) is the coarseness of the
%       i'th variable as determined by the variable ranges and
%       BITS(i).  ROUND contains the absolute indices of the
%       X which where rounded due to finite BIT length.
%
%	Copyright (c) 1993 by the MathWorks, Inc.
%	Andrew Potvin 1-10-93.

% Remark: what about handling case where length(bits)~=length(vlb)?


lchrom = sum(bits);
coarse = (vub-vlb)./((2.^bits)-1);
[x_row,x_col] = size(x);

gen = [];
if ~isempty(x),
   temp = (x-ones(x_row,1)*vlb)./ ...
          (ones(x_row,1)*coarse);
   b10 = round(temp);
   % Since temp and b10 should contain integers 1e-4 is close enough
   nround = find(b10-temp>1e-4);
   gen = b10to2(b10,bits);
end

% end encode


function [x,coarse] = decode(gen,vlb,vub,bits)
%DECODE Converts from binary Gray code to variable representation.
%	[X,COARSE] = DECODE(GEN,VLB,VUB,BITS) converts the binary 
%       population GEN to variable representation.  Each individual 
%       of GEN should have SUM(BITS).  Each individual binary string
%       encodes LENGTH(VLB)=LENGTH(VUB)=LENGTH(BITS) variables.
%       COARSE is the coarseness of the binary mapping and is also
%       of length LENGTH(VUB).
%
%  this *.m file created by combining "decode.m" from the MathWorks, Inc.
%  originally created by Andrew Potvin in 1993, with "GDECODE.FOR" written 
%  by William A. Crossley in 1996.
%	
%	William A. Crossley, Assoc. Prof. School of Aero. & Astro.
%  Purdue University, 2001
%
%  gen is an array [population size , string length], each row is one individual's chromosome
%  vlb is a row vector [number of parameters], each entry is the lower bound for a variable
%  vub is a row vector [number of parameters], each entry is the upper bound for a variable
%  bits is a row vector [number of parameters], each entry is number of bits used for a variable
%  

no_para = length(bits); % extract number of parameters using number of rows in bits vector
npop = size(gen,1);		% extract population size using number of rows in gen array
x = zeros(npop, no_para);  % sets up x as an array [population size, number of parameters]
coarse = zeros(1,no_para); % sets up coarse as a row vector [number of parameters]

for J = 1:no_para,  % extract the resolution of the parameters
	coarse(J) = (vub(J)-vlb(J))/(2^bits(J)-1);	% resolution of parameter J
end

for K = 1:npop, 	% outer loop through each individual (there may be a more efficient way to operate on the
                  % gen array) BC 10/10/01
	sbit = 1;		% initialize starting bit location for a parameter
	ebit = 0;		% initialize ending bit location
   
   for J = 1:no_para,	% loop through each parameter in the problem
   	ebit = bits(J) + ebit;	% pick the end bit for parameter J
		accum = 0.0;				% initialize the running sum for parameter J
      ADD = 1;						% add / subtract flag for Gray code; add if(ADD), subtract otherwise
      for I = sbit:ebit,			% loop through each bit in parameter J
         pbit = I + 1 - sbit;		% pbit determines value to be added or subtracted for Gray code
         if (gen(K,I))					% if "1" is at current location
            if (ADD)						% add if appropriate
               accum = accum + (2.0^(bits(J)-pbit+1) - 1.0);
               ADD = 0;					% next time subtract
            else
               accum = accum - (2.0^(bits(J)-pbit+1) - 1.0);
               ADD = 1;					% next time add
            end
         end
      end								% end of I loop through each bit
      x(K,J) = accum * coarse(J) + vlb(J);			% decoded parameter J for individual K
      sbit = ebit + 1;										% next parameter starting bit location
   end						% end of J loop through each parameter
end 					% end of K loop through each individual

%end gdecode


function [new_gen,mutated] = mutate(old_gen,Pm)
%MUTATE Changes a gene of the OLD_GEN with probability Pm.
%	[NEW_GEN,MUTATED] = MUTATE(OLD_GEN,Pm) performs random
%       mutation on the population OLD_POP.  Each gene of each
%       individual of the population can mutate independently
%       with probability Pm.  Genes are assumed possess boolean
%       alleles.  MUTATED contains the indices of the mutated genes.
%
%	Copyright (c) 1993 by the MathWorks, Inc.
%	Andrew Potvin 1-10-93.

mutated = find(rand(size(old_gen))<Pm);
new_gen = old_gen;
new_gen(mutated) = 1-old_gen(mutated);

% end mutate


function [new_gen,nselected] = tourney(old_gen,fitness)
%TOURNEY Creates NEW_GEN from OLD_GEN, based on tournament selection.
%	 [NEW_GEN,NSELECTED] = TOURNEY(OLD_GEN,FITNESS) selects
%        individuals from OLD_GEN by competing consecutive individuals
%	 after random shuffling.  NEW_GEN will have the same number of
%	 individuals as OLD_GEN.
%        NSELECTED contains the number of copies of each individual
%	 that survived.  This vector corresponds to the original order
%	 of OLD_GEN.
%
%	 Created on 1/21/96 by Joel Grasmeyer

% Initialize nselected vector and indices of old_gen
new_gen = [];
nselected = zeros(size(old_gen,1),1);
i_old_gen = 1:size(old_gen,1);

% Perform two "tournaments" to generate size(old_gen,1) new individuals
for j = 1:2,

  % Shuffle the old generation and the corresponding fitness values
  [old_gen,i_shuffled] = shuffle(old_gen);
  fitness = fitness(i_shuffled);
  i_old_gen = i_old_gen(i_shuffled);

  % Keep the best of each pair of individuals
  index = 1:2:(size(old_gen,1)-1);
  [min_fit,i_min] = min([fitness(index);fitness(index+1)]);
  selected = i_min + [0:2:size(old_gen,1)-2];
  new_gen = [new_gen; old_gen(selected,:)];

  % Increment counters in nselected for each individual that survived
  temp = zeros(size(old_gen,1),1);
  temp(i_old_gen(selected)) = ones(length(selected),1);
  nselected = nselected + temp;

end

% end tourney


function [new_gen,index] = shuffle(old_gen)
%SHUFFLE Randomly reorders OLD_GEN into NEW_GEN.
%	 [NEW_GEN,INDEX] = MATE(OLD_GEN) performs random reordering
%        on the indices of OLD_GEN to create NEW_GEN.
%	 INDEX is a vector containing the shuffled row indices of OLD_GEN.
%
%	 Created on 1/21/96 by Joel Grasmeyer

[junk,index] = sort(rand(size(old_gen,1),1));
new_gen = old_gen(index,:);

% end shuffle


function [new_gen,sites] = uniformx(old_gen,Pc)
%UNIFORMX Creates a NEW_GEN from OLD_GEN using uniform crossover.
%	  [NEW_GEN,SITES] = UNIFORMX(OLD_GEN,Pc) performs uniform crossover
%         on consecutive pairs of OLD_GEN with probability Pc.
%	  SITES shows which bits experienced crossover.  1 indicates
%	  allele exchange, 0 indicates no allele exchange.  SITES has
%	  size(old_gen,1)/2 rows.
%
%  	  Created 1/20/96 by Joel Grasmeyer

new_gen = old_gen;
sites = rand(size(old_gen,1)/2,size(old_gen,2)) < Pc;
for i = 1:size(sites,1),
  new_gen([2*i-1 2*i],find(sites(i,:))) = old_gen([2*i 
2*i-1],find(sites(i,:)));
end

% end uniformx


