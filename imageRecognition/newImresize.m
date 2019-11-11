function tst = newImresize(M)

% M = single(imread('tst.png')); % this  part's the same
[H,W,~]=size(M); % sometimes it's easier to use "H=height" and "W=width" to keep things in order

% generate X and Y grids for old and new sizes
x_orig=linspace(0,1,W);
y_orig=linspace(0,1,H); % get X and Y grids for original pixels
[X_orig,Y_orig]=meshgrid(x_orig,y_orig);
scale=1024/max(H,W);
W_new=round(W*scale); H_new=round(H*scale);
x_new=linspace(0,1,round(W_new));
y_new=linspace(0,1,round(H_new));
[X_new,Y_new]=meshgrid(x_new,y_new);

tst=zeros(H_new,W_new,2,'uint8');
tst(:,:,1) = uint8(interp2(X_orig,Y_orig,M,X_new,Y_new,'cubic'));
tst(:,:,2) = unit8(imresize(M,[1024 1024]));

end