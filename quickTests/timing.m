function [fList_Output,FList_Output] = timing(this,z)

[X,Y] = meshgrid(this.grid_fxList,this.grid_fyList);

% Copied function to mimic the Heaviside.m function from SymMath Toolbox without using it
hvsd = @(x)(0.5*(x == 0) + (x > 0));

% % Phase aquired by traveling a distance z forward in time
% % phaseVec = complex(zeros(this.grid_npts));
% phaseVec = exp( 1i * ( (2 * pi) / (this.grid_lambda) ) *...
%     sqrt(1 - X.^2 * this.grid_lambda^2 - Y.^2 * this.grid_lambda^2) * z )...
%     .* hvsd(1 - (X.^2 * this.grid_lambda^2 + Y.^2 * this.grid_lambda^2) );
%     % Heaviside is to keep waves from becoming even evanescent
%
% % This first transform takes us into the spatial frequencies where we apply
% % the propagation phase
% %%% For consistancy I have not obeyed camelCase here %%%
% FList_Output =...
%     this.grid_dx .* this.grid_dy .*...
%     ifftshift(fft2(fftshift(this.field_fList)))...
%     .* phaseVec;
%
% % This second transform takes us back to the real space
% fList_Output =...
%     (1 / (this.grid_dx .* this.grid_dy) ).*...
%     ifftshift(ifft2(fftshift(FList_Output)))...
%     ;


%%



A = ( (2 * pi) / (this.grid_lambda) ) * sqrt(1 - X.^2 * this.grid_lambda^2 - Y.^2 * this.grid_lambda^2) * z;

tic
for ii = 1:10
    B = conj(exp( 1i * A ));
end
toc/10

tic
for ii = 1:10
    C = exp( -1i * A );
end
toc/10


%%
FList_Output = [];
fList_Output = [];


end