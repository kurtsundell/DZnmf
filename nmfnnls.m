function [source_PDP,Weightings,numIter,tElapsed,finalResidual]=nmfnnls(Basin_PDP,nsources,option)
%Non-negative Matrix Factorization based on Non-negative Least Squares optimization
%X=AY, s.t. X,A,Y>=0.
%
%Definition:
%     [A,Y,numIter,tElapsed,finalResidual]=nmfnnls(Basin_PDP,k)
%     [A,Y,numIter,tElapsed,finalResidual]=nmfnnls(Basin_PDP,k,option)
%Inputs:
% Basin_PDP: non-negative matrix of detrital PDPs  to factorize, each column is a 
%   sample, and each row is a relative probability for a given age.
% nsources: scalar, number of sources (factors) to factorize Basin_PDP into.
% option: struct:
%   option.iterations: max number of interations. The default is 10000.
%   option.dis: boolen scalar, It could be 
%        false: not display information ,
%        true: display (default).
%   option.residual: the threshold of the fitting residual to terminate. 
%        If the ||X-XfitThis||<=option.residual, then halt. The default is 1e-8.
%   option.tof: if ||XfitPrevious-XfitThis||<=option.tof, then halt. The default is 2e-16.
%
%OUTPUTS
% source_PDP: the basis matrix (factorized sources).
% Weightings: the coefficient matrix (factorized weighting fuctions).
% numIter: scalar, the number of iterations.
% tElapsed: scalar, the computing time used.
% finalResidual: scalar, the fitting residual.
% References:
%  [1]\bibitem{NMF_ANLS_Kim2008}
%     H. Kim and H. Park,
%     ``Nonnegative matrix factorization based on alternating nonnegativity constrained least squares and active set method,''
%     {\it SIAM J. on Matrix Analysis and Applications},
%     vol. 30, no. 2, pp. 713-730, 2008.
%  [2]\bibitem{NMF_Sparse_Kim2007}
%     H. Kim and H. Park,
%     ``Sparse non-negatice matrix factorization via alternating non-negative-constrained least squares for microarray data analysis,''
%     {\it Bioinformatics},
%     vol. 23, no. 12, pp. 1495-1502, 2007.
tStart=tic;
%Set stopping criteria
optionDefault.iterations=10000;
optionDefault.display=true;
optionDefault.residual=1e-8;
optionDefault.tof=2e-16;
if nargin<3
   option=optionDefault;
else
    option=mergeOption(option,optionDefault);
end


%initialize variables
[r,c]=size(Basin_PDP); % c is # of samples, r is # of features
Weightings=rand(nsources,c); %initiallize Weightings matrix
XfitPrevious=Inf; %initiallize fit 
source_PDP=rand(r,nsources); %initiallize source_PDP matrix

%Initialize waitbar
if option.display
     h=waitbar(0,'Calculating','CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
     setappdata(h,'canceling',0)
end
%repeat NMF iterations until the specified number of iterations is reached
...or until another stopping criterion is reached
for i=1:option.iterations 
    source_PDP=kfcnnls(Weightings',Basin_PDP'); 
    source_PDP=source_PDP';
    source_PDP=normc(source_PDP);
    
    %normalize sum of source_PDP to 1
    sumA=cumsum(source_PDP);
    max_A=sumA(r,:);
    source_PDP=bsxfun(@rdivide,source_PDP,max_A);
    
    Weightings=kfcnnls(source_PDP,Basin_PDP);
    %normalize sum of Weightings to 1
    sumY=cumsum(Weightings);
    max_Y=sumY(nsources,:);
    Weightings=bsxfun(@rdivide,Weightings,max_Y);
    
   if mod(i,20)==0 || i==option.iterations
        
        XfitThis=source_PDP*Weightings;
        iterationResidual=norm(XfitPrevious-XfitThis, 'fro')%calculate change in 
                            ...factorized source_PDP for this iteration
        XfitPrevious=XfitThis; %update X reconstruction
        currentResidual=norm(Basin_PDP-XfitThis,'fro'); %calculate final residual
                            ... for this iteration
                                
                        
        % End iterations if one of the ending criterion are met
        if option.tof>=iterationResidual || option.residual>=currentResidual || i==option.iterations
            s=sprintf('NNLS based NMF successes! \n # of iterations is %0.0d. \n The final residual is %0.4d.',i,currentResidual);
            disp(s);
            numIter=i;
            finalResidual=currentResidual;
            break;
        end
   end
   if option.display %break iterations and delete waitbar if cancelled
       if getappdata(h,'canceling')
           delete(h)
           break
       end
       waitbar(i/option.iterations)%update waitbar
   end
   
end
if option.display%delete waitbar after calculations are complete
    delete(h)
end
tElapsed=toc(tStart);
end

