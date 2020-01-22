function [PDP,Weightings,finalResidual, coefficient_count,R2, numIter]=DZnmf_loop(Basin_PDP, nsources)
%Calculates the source characteristics for all NMF ranks from 2 to nsources
%INPUT VARIABLES:
%Basin_PDP: a-by-n matrix of n sink mixture distributions (arranged in columms)
%   of a relative probabilities measured at regular intervals over the age
%   interval [x, y]
%nsources: number of sources to reconstruct
%OUTPUT VARIABLES
%PDP: a-by-m-by-(nsources-1) matrix of m source mixture distributions (arranged in columns)
%   of a relative probabilties measured at regular intervals over the age
%   range [z,y]
%Weightings: m-by-n-by-(nsources-1) matrix of weighting functions describing the mixture of
%   PDP sources that reconstruct Basin_PDP sink samples
%finalResidual: the final misfit between reconstructed and input Basin_PDP
%   matrixes
%coefficient_count:NMF ranks tried
%R2: 1-by-n-by-(nsources-1) matrix of Cross-correlation coefficients between
%   reconstructed and input Basin_PDP matrices
%numIter: number of iterations required to reach ending conditions

global iter_num
global cancel

%Initialize variables
rows=size(Basin_PDP, 1);
columns=size(Basin_PDP,2);
PDP=zeros(rows,2, nsources-1); %initialize output 3D array 
Weightings=zeros(2,columns, nsources-1); %initialize output wieghting 3D matrix 
coefficient_count=2;

%Run NMF
for i=2:nsources
	if cancel == 1
		err_dlg=errordlg('Calculation stopped','');
		waitfor(err_dlg);
		break
	else
		iter_num = i;
		[PDP(:,:,i-1),Weightings(:,:,i-1),numIter(i-1),tElapsed,finalResidual(i-1)]=DZnmf(Basin_PDP,i);
		PDP(:,i+1,i-1)=zeros(rows,1);
		Weightings(i+1,:,i-1)=zeros(1,columns);
		coefficient_count(i)=coefficient_count(i-1)+1;
	end
end

if cancel == 0	

%Trim Matricies
PDP=PDP(:,1:nsources,:);
Weightings=Weightings(1:nsources,:,:);
coefficient_count=coefficient_count(:,1:nsources-1);

%Compare output PDPs
for i=1:nsources-1
    PDP_output(:,:,i)=PDP(:,:,i)*Weightings(:,:,i);
end
CDF = cumsum(PDP_output);
B=Basin_PDP-mean(Basin_PDP);
trash3=sum(B.*B);
trash33SQRT=sqrt(trash3);
for i=1:nsources-1
    A=PDP_output(:,:,i) - mean(PDP_output(:,:,i));
    trash1=sum(A.*B);
    trash2=sum(A.*A);
    trash22SQRT=sqrt(trash2);
    R2(:,:,i) = trash1./(trash22SQRT.*trash33SQRT);
    %R2(:,:,i) = (sum(trash1))/(sqrt((sum(trash2))*(sum(trash3))))^2;
end

for i=1:nsources-1
    R2_mean(i)=mean(R2(1,:,i));
    R2_std(i) = std(R2(1,:,i));    
end

end
    
