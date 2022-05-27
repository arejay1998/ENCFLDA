function [ pearsSimilarity ] = pearsSim( vectorA, vectorB )
    pearsSimilarityMatrix = 0.5 + 0.5 * corrcoef(vectorA, vectorB);
    pearsSimilarity = pearsSimilarityMatrix(1,2);
end
