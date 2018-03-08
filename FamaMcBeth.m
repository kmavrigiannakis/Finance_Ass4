clear all;
clc;

load('size.mat');
sizeIndustry= Size;
clear Size;

load('booktomarket.mat');

load('DATAFinanceExercise2.mat');
monthlyReturns = AverageValueWeightedReturnsMonthly;
monthlyReturns(isnan(monthlyReturns)==1 )=0;
clear AverageValueWeightedReturnsMonthly;

load('STR.mat');
load('momentumMatrix.mat');

load('DATEBOOKTOMARKET.mat');

%Creating monthly book to market matrix
bookToMarketMonthly = zeros( 12*size(booktomarket,1), size(monthlyReturns,2 ) );
temp=1;
positionInBookToMarketMatrix =1;
for i=1:12*size(booktomarket,1)
    if(temp==13)
        temp=1;
        positionInBookToMarketMatrix =positionInBookToMarketMatrix +1;
    end
    bookToMarketMonthly(i,:)=booktomarket(positionInBookToMarketMatrix,:) ;
    temp=temp+1;
end

%Breaking Matrices so as that all to have the same dimentions
sizeIndustry = sizeIndustry(15:end ,:);
bookToMarketMonthly = bookToMarketMonthly(20:end,:);
monthlyReturns = monthlyReturns(14:end,:);

%Regression 
numMonths = size(momentum,2);
numAssets = size(monthlyReturns,2);

%Reshaping appropriately the matrices for the regression
momentum = reshape(momentum, numMonths,1 );
STR = reshape(STR, numMonths,1 );

%Initialization of beta matrix for each asset
betaMatrix = zeros(numAssets,4 );
%rsquareMatrix = zeros( )

%1st Step of Fama MacBeth
%Regression for each asset of its return on the 4 factors
for i=1:numAssets
 %i=1;
    %Creating matrix of regressors
    regressors = [ STR, momentum , sizeIndustry( :,i ) , bookToMarketMonthly(:,i) ];
 
    regStats = regstats( monthlyReturns(:,i),regressors ,'linear');
    betaMatrix(i,:) = regStats.beta(2:end,:);
    %regStats.rsquare  
    
end

%Estimation of average returns for each asset
meanReturns = reshape(  mean( monthlyReturns,1) , numAssets,1  );

% %2nd Step 
% %Regression of average return of each asset, on coefficients estimated
% %on first step
%     regStats = regstats(meanReturns, betaMatrix,'linear');
%     factorPremia= regStats.beta; 

monthlyReturns2 = reshape(monthlyReturns, numAssets, numMonths);
%Initialization of matrices of coefficients and statistics
factorPremia= zeros(numMonths,5);
rsquareVector = zeros(numMonths,1);
tstatistic = factorPremia;
pvalue= factorPremia;
matrixNwse=zeros(numMonths,5);

for t=1:numMonths
%for t=1:45 
    %for i=1:numAssets 
   
    regStats2 = regstats(monthlyReturns2(:,t), betaMatrix,'linear'); 
       % regStats2 = regstats(monthlyReturns(t,:), a,'linear');      
    factorPremia(t,:) = regStats2.beta; 
    rsquareVector(t) = regStats2.rsquare;
    tstatistic(t,:) = regStats2.tstat.se;
    pvalue(t,:)= regStats2.tstat.pval;
    matrixNwse(t,:) =nwse(regStats2.r , betaMatrix  );
    
    %end
end

averageFactorPremia = nanmean(factorPremia);
averageRsquare = nanmean(rsquareVector);


