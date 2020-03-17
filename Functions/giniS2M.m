load('./Data/CriteriaAnalysis/AsymmetricParameterSet/sol5_3_1')

countgini = 1;
for i = 1:5
    genecount = sol{35}.samp(i,:);
    G(countgini) = gini(ones(1000,1),genecount);
    countgini = countgini+1;
end

load('./Data/CriteriaAnalysis/5nodes/sol10005_1_55')
param = 4;

countgini = 1;
for i = 1:5
    genecount = sol{param}.samp(i,:);
    GII(countgini) = gini(ones(1000,1),genecount);
    countgini = countgini+1;
end
