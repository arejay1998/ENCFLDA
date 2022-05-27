function [ cosSimilarity ] = cosSim1( vectorA, vectorB )
    %ע��vectorA��vectorB����������
    num = vectorA * vectorB';
    denom = norm(vectorA) * norm(vectorB);
    cosSimilarity = 0.5 + 0.5 * (num./denom);
end
