function timing(scriptName)

N = 2^10;
M = 51;

largeMat = zeros(N,N);
smallMat = ones(M,M);

tic
for ii = 1:2625
largeMat(...
    (size(largeMat,1)/2-size(smallMat,1)/2):(size(largeMat,1)/2+size(smallMat,1)/2)-1,...
    (size(largeMat,1)/2-size(smallMat,1)/2):(size(largeMat,1)/2+size(smallMat,1)/2)-1 ) =...
    largeMat(...
    (size(largeMat,1)/2-size(smallMat,1)/2):(size(largeMat,1)/2+size(smallMat,1)/2)-1,...
    (size(largeMat,1)/2-size(smallMat,1)/2):(size(largeMat,1)/2+size(smallMat,1)/2)-1 ) + smallMat;
end
toc

end