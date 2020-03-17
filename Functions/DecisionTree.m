%function generatign decision tree, where rare coordinated high parameter
%sets are class 1, the other parameter sets in class 0 
%
%INPUT:
%
%I:         rare coordinated high parameter sets
%R_an:      matrix of all independent parameters - exlcding the dependent
%           parameter k
%
%OUTPUT:
%
%decision tree plot

function DecisionTree(I,R_an)

oversamp = round((length(R_an)-length(I))/length(I))-1;

X = R_an;

Y = zeros(length(R_an),1);
Y(I) = 1;

tree = fitctree(X,Y);
view(tree,'Mode','graph');

end