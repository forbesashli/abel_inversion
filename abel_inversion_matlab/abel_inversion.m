function [ f_rec , X ] = abel_inversion(h,R,upf,plot_results,lsq_solve)
% This function calculates a Fourier-based Abel inversion based on the
% method described in [1]. Assuming cylindrical symmetry, the distribution 
% function f_rec can be reconstructed from the measured profile h .
%
% The main idea is a Fourier-series-like expansion of the unknown
% distribution function f(r) into 
% f(r) = sum_{n=lof}^{upf} (A_n * f_n(r))                               (1)
% where the lower frequency is set to 1 and the upper frequency upf is
% important for noise-filtering.
%
%
% INPUT:  H - dataset to be analyzed (1xN - matrix). 
%           For optimum results, H(1) should be the center and H(N) the 
%           edge of the investigated object and H(N) should be approx zero.
%           If no input H is given, a sampla dataset will be created.
%         R - radius of given system (i.e. distance between H(1) and H(N))
%         UPF  - upper frequency limit. Defines the number of cosinus
%           expansions (and critically determines computation time). 
%           Choosing a very low value (e.g. 4) results in a low-pass 
%           filtering effect, reducing noise (but also potential features) 
%         PLOT_RESULTS - set to 1 to plot the results
%         LSQ_SOLVE - set to 1 to use LSQCURVEFIT instead of algebraically
%           solving the least square problem - maybe useful if the inversion
%           of a singular matrix ocurrs (I am not sure if this might happen)
%
% OUTPUT: F_REC - reconstructed density profile
%         X - x-vector containing spatial coordinates for F_REC
%
%    
%     [1] G. Pretzler, Z. Naturforsch. 46a, 639 (1991)
%
% See also: COMPUTE_EXPANSION, GENERATE_TEST_DATA, SOLVE_LSQ
%
%                                         written by C. Killer, Sept. 2013


%% format data / generate sample data

% create sample data based on a polynomial distribution function
% if no data input is given
if ~exist('h', 'var') || isempty(h)
    [X,h,R]=generate_test_data;
    plot_results=1;
else
    X=linspace(0,R,length(h))';
end 

% default value for number of expansion elements
if ~exist('upf', 'var'); upf=10; end; 

% avoiding problems if flags are not given in input
if ~exist('plot_results', 'var'); plot_results=0; end; 
if ~exist('lsq_solve', 'var'); lsq_solve=0; end; 

%% calculate series expansions fn and corresponding integrals hn

[fn,hn] = compute_expansion( X,upf,R );

%% solve equation system A*L=B for the amplitudes A
if lsq_solve ~= 1
    
    B = zeros(1,upf+1); L = zeros(upf+1,upf+1);  %create arrays
    
    for k=1:upf+1

        for l=1:upf+1
            L(l,k)=2.*sum(hn(:,k).*hn(:,l)); % January 2016: added factor 2
        end  

        B(k)=sum(h(:).*hn(:,k));
    end
    
    A=B/L;
    
else
       
    x0=1*ones(upf+1,1);       % guess some initial values for optimisation     
    A=solve_lsq(h,hn,x0);     % solve for amplitudes A
end

%% final stage: calculate the resulting distribution profile

% create vector for resulting reconstructed distribution
f_rec=zeros(length(h),1);       

% special case for n=0 (where f_0(r) = 1)
f_rec = f_rec + A(1)*1;

% iterate eq. (1) for n=1:upf
for c=2:upf+1
    f_rec = f_rec + A(c).*fn(:,c);
end

if plot_results==1    
    figure; % normalized profiles for better comparison
    set(gca,'linewidth',1.5,'fontsize',16)
    hold on; 
    plot(X,h./max(h),'b','Linewidth',1.5); 
    plot(X,f_rec./max(f_rec),'k','Linewidth',1.5); 
    grid on; box on; 
    title(sprintf('number of cos-expansions: %i',upf))
    legend('measured profile','reconstructed distribution','Location','SouthWest')
end

