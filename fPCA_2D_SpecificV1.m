function [PrimaryComponent, PCDirection, PCVariance]=...
    fPCA_2D_SpecificV1(Imagery2D,Para, VarOrNum, Style)
%VarOrNum���ж�Para���������ɷַ��������ȣ��������ɷָ���
%                   VarOrNum=1 ParaΪ���ɷַ��������ȣ�VarOrNum=0 ParaΪ���ɷָ���
%Style 0 normal PCA
%        1 regularized PCA
%

[SamNum,Dim] = size(Imagery2D);%ʹ�õ�Դ����

CoVar = cov(Imagery2D);
[V,D] = eig(CoVar) ;


% if rank(CoVar) == Dim
% [V,D] = eig(CoVar) ;
% else
%     [V,D] = eig(CoVar+0.00001*eye(Dim,Dim)) ; % to avoid the singularity of matrix
% end
    
for i=1:Dim
    EigenValueTemp(i)=D(i,i);%������ֵ��ȡ��һ������
end
[EigenValue,Index]= sort(EigenValueTemp,'descend');%����ֵ��С��������
EigenMatrix=V(:,Index);

if  VarOrNum==1
     EigenValueSum=sum(EigenValue);
     TempEigenValueSum = 0;
     EigenVectorNum = 0;
     for i=1:Dim
            Temp=EigenValue(i);
            TempEigenValueSum = TempEigenValueSum+Temp;
             EigenValueWeigeht=TempEigenValueSum/EigenValueSum;
             if EigenValueWeigeht>Para
                break;    
             end
             EigenVectorNum=EigenVectorNum+1;             
     end
 end

 if  VarOrNum==0
     EigenVectorNum=Para;
 end

 if Style == 1
     for i = 1:size(EigenMatrix,2)
     TempNorm(1,i) = norm(EigenMatrix(:,i));
     end
     EigenMatrix = EigenMatrix./repmat(TempNorm,size(EigenMatrix,1),1);
 end

PCDirection = EigenMatrix(:,1:EigenVectorNum);
PCVariance = EigenValue(EigenVectorNum);
PrimaryComponent = Imagery2D*PCDirection;
%TotalPriComp = Imagery2D*EigenMatrix;
