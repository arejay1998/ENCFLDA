function [ cosSimilarity ] = cosSim1( vectorA, vectorB )
    %注意vectorA和vectorB都是行向量
    num = vectorA * vectorB';
    denom = norm(vectorA) * norm(vectorB);
    cosSimilarity = 0.5 + 0.5 * (num./denom);
end
