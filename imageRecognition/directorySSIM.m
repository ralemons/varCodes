clear


for ii = 1:5
    fileName = ['C:\Users\rlemons\Documents\Pictures\Test Roses\r',num2str(ii),'.jpg'];
    img{ii} = imresize(imread(fileName),[300 300]);
    if size(img{ii},3) == 3
        img{ii} = rgb2gray(img{ii});
    end
end

tst = repmat(img,1,10);

N = length(tst);

tmp = zeros(N);

for zz = 1:5
tic
for jj = 1:N
    for ii = jj+1:N
        tmp(jj,ii) = slimSSIM(tst{jj},tst{ii});
    end
end
toc
end