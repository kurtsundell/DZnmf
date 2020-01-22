function [source_PDP,Weightings,numIter,tElapsed,finalResidual]=nmfrule(X,nsources,option)
% NMF based on multiple update rules: X=AY, s.t. X,A,Y>=0.
% Definition:
%     [A,Y,numIter,tElapsed,finalResidual]=nmfrule(X,k)
%     [A,Y,numIter,tElapsed,finalResidual]=nmfrule(X,k,option)
% X: non-negative matrix, dataset to factorize, each column is a sample, and each row is a feature.
% k: number of clusters.
% option: struct:
% option.distance: distance used in the objective function. It could be
%    'ls': the Euclidean distance (defalut),
%    'kl': KL divergence.
% option.iter: max number of interations. The default is 1000.
% option.dis: boolen scalar, It could be 
%     false: not display information,
%     true: display (default).
% option.residual: the threshold of the fitting residual to terminate. 
%    If the ||X-XfitThis||<=option.residual, then halt. The default is 1e-4.
% option.tof: if ||XfitPrevious-XfitThis||<=option.tof, then halt. The default is 1e-4.
% A: matrix, the basis matrix.
% Y: matrix, the coefficient matrix.
% numIter: scalar, the number of iterations.
% tElapsed: scalar, the computing time used.
% finalResidual: scalar, the fitting residual.
% References:
%  [1]\bibitem{NMF_Lee1999}
%     D.D. Lee and S. Seung,
%     ``Learning the parts of objects by non-negative matrix factorization,''
%     {\it Science},
%     vol. 401, pp. 788-791, 1999.
%  [2]\bibitem{NMF_Lee2001}
%     D.D. Lee and S. Seung,
%     ``Algorithms for non-negative matrix factorization,''
%     {\it Advances in Neural Information Processing Systems}, 
%     vol. 13, pp. 556-562, 2001.
%%%%
% Copyright (C) <2012>  <Yifeng Li>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% Contact Information:
% Yifeng Li
% University of Windsor
% li11112c@uwindsor.ca; yifeng.li.cn@gmail.com
% May 01, 2011
%%%%

tStart=tic;
optionDefault.distance='ls';
optionDefault.iter=100000;
optionDefault.dis=true;
optionDefault.residual=1e-8;
optionDefault.tof=2e-16;
if nargin<3
   option=optionDefault;
else
    option=mergeOption(option,optionDefault);
end

% iter: number of iterations
[r,c]=size(X); % c is # of samples, r is # of features
Weightings=rand(nsources,c);
% Y(Y<eps)=0;
Weightings=max(Weightings,eps);

%normalize sum of Weightings to 1
sumY=cumsum(Weightings);
max_Y=sumY(nsources,:);
Weightings=bsxfun(@rdivide,Weightings,max_Y);

source_PDP=X/Weightings;
% A(A<eps)=0;
source_PDP=max(source_PDP,eps);

%normalize sum of source_PDP to 1
sumA=cumsum(source_PDP);
max_A=sumA(r,:);
source_PDP=bsxfun(@rdivide,source_PDP,max_A);
XfitPrevious=Inf;
for i=1:option.iter
    switch option.distance
        case 'ls'
            source_PDP=source_PDP.*((X*Weightings')./(source_PDP*(Weightings*Weightings')));
%             A(A<eps)=0;
                source_PDP=max(source_PDP,eps);
            Weightings=Weightings.*((source_PDP'*X)./(source_PDP'*source_PDP*Weightings));
%             Y(Y<eps)=0;
                Weightings=max(Weightings,eps);
        case 'kl'
            source_PDP=source_PDP.*(((X./(source_PDP*Weightings))*Weightings')./(ones(r,1)*sum(Weightings,2)'));
            source_PDP=max(source_PDP,eps);
            Weightings=Weightings.*((source_PDP'*(X./(source_PDP*Weightings)))./(sum(source_PDP,1)'*ones(1,c)));
            Weightings=max(Weightings,eps);
        otherwise
            error('Please select the correct distance: option.distance=''ls''; or option.distance=''kl'';');
    end
    if mod(i,100)==0 || i==option.iter
        if option.dis
            disp(['Iterating >>>>>> ', num2str(i),'th']);
        end
        XfitThis=source_PDP*Weightings;
        fitRes=matrixNorm(XfitPrevious-XfitThis);
        XfitPrevious=XfitThis;
        curRes=norm(X-XfitThis,'fro');
        if option.tof>=fitRes || option.residual>=curRes || i==option.iter
            s=sprintf('Mutiple update rules based NMF successes! \n # of iterations is %0.0d. \n The final residual is %0.4d.',i,curRes);
            disp(s);
            numIter=i;
            finalResidual=curRes;
            break;
        end
    end
end
tElapsed=toc(tStart);
end
