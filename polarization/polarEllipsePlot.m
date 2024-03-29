function [varargout] = polarEllipsePlot(S,fullImage,div,arrowScale,arrowPos,fignum,chooseROI,macroPixSize,logMat)
% This function makes use of the Image Processing Toolbox from MATLAB. If
% you do not have this toolbox then you cannot use this function unless you
% want to rebuild the functionality of 'drawcircle' and 'inROI' yourself.

% X and Y locations of the arrows
X = 1:div:size(S,2);
Y = 1:div:size(S,1);

% % Scales image and places it on grid for the ellipses to plot on
f = figure(fignum);
clf
% imagesc([0 size(S,2)+1],[0 size(S,1)+1], fullImage );
xlim([0 size(S,2)+1]);
ylim([0 size(S,1)+1]);
axis square off

% ROI choice and plotting
if chooseROI
    h = drawrectangle('DrawingArea','unlimited');
    customWait(h);
    h.Visible = 'off';
else
    h = images.roi.Rectangle(gca,...
        'Position',[0,0,size(S,2)+1,size(S,1)+1],...
        'Visible','off');
end


% Takes X and Y into a list of ordered points
XY(:,1) = repmat(X,1,numel(X));
XY(:,2) = sort(repmat(Y,1,numel(Y)));

% Built in function to find which points are in the ROI
XYinROI = inROI(h,XY(:,1),XY(:,2));
XY = XY(XYinROI,:);

% if ~exist('macroPixSize','var') && ~exist('logMat','var')
if ~exist('logMat','var')
    logMat = ones(size(fullImage));
end

% Sets size of figure

hold on

% % Holds info on whether arrow flips direction based on e>1
% flip = 0;

counter = 1;
co = [255, 255, 255]/255;


% Loop through all the ellipses that need to be created
for ii = 1:size(XY,1)
    
    if logMat((XY(ii,2)*macroPixSize.V),(XY(ii,1)*macroPixSize.H)) == 1
        
        
        % Find stokes parameters
        s0 = S(XY(ii,2),XY(ii,1),1);
        s1 = S(XY(ii,2),XY(ii,1),2);
        s2 = S(XY(ii,2),XY(ii,1),3);
        s3 = S(XY(ii,2),XY(ii,1),4);
        
        % Generate ellipse values based on stokes parameters
        e = sqrt( (2 * sqrt(s1^2+s2^2)) / (1 + sqrt(s1^2+s2^2)) );
        theta = pi-atan2(s2,s1)/2;
        
        % Major and minor axis of the ellipse
        a = (s0/2)*div;
        b = a.*sqrt(1-e.^2);
        
        % If e>1 then b goes pure imaginary. This sets b as the negative value
        % (the ellipse is pointing the other way) and changes the direction of
        % the arrow to maintain correct handedness
        if ~isreal(b)
            b = -imag(b);
            s3 = -s3;
            flip = 1;
        else
            flip = 0;
        end
        
        % Generate the ellipse points centered on (0,0)
        nPts = 50;
        A = a .* cos(0:pi/nPts:2*pi);
        B = b .* sin(0:pi/nPts:2*pi);
        
        
        % Rotate the ellipse around (0,0)
        rotMat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        xy = [A;B]' * rotMat;
        A = xy(:,1);
        B = xy(:,2);
        
        % Place ellipse in correct place ie. (0,0) -> (X,Y)
        A = A + XY(ii,1);
        B = B + XY(ii,2);
        
        % Where to plot arrow around the ellipse
        if flip == 0
            [~,top] = min( abs((0:pi/nPts:2*pi) - pi/2) );
        else
            [~,top] = min( abs((0:pi/nPts:2*pi) - 3*pi/2) );
        end
                
        % Which direction to plot the arrow
        if sign(s3) == 1 || sign(s3) == 0
            xArrowHead = [A(top-1) A(top)];
            yArrowHead = [B(top-1) B(top)];
%             co = [255, 255, 255]/255;
        else
            xArrowHead = [A(top) A(top-1)];
            yArrowHead = [B(top) B(top-1)];
%             co = [255, 73, 41]/255;
        end
        
        % Keeps the warning from showing up about imaginary numbers
        warning(''); %#ok<WLAST>
        
        % Colors for plotting
        %     co = [0.8500 0.3250 0.0980];
        %     coWarn = [0.7500 0 0.7500];
        
        if a > 0.2
            
            p(counter) = plot(A,B,'Color',co,'linewidth',1); % Ellipse
            if isempty(arrowPos) && isempty(arrowScale * s0)
                ar(counter) = arrowh(xArrowHead,yArrowHead,co); % Arrow
            elseif isempty(arrowPos) && ~isempty(arrowScale * s0)
                ar(counter) = arrowh(xArrowHead,yArrowHead,co,arrowScale * s0); % Arrow
            else
                ar(counter) = arrowh(xArrowHead,yArrowHead,co,arrowScale * s0,arrowPos); % Arrow
            end
            
        end
        
%         if flip == 1
%             p(counter).Color = coWarn;
%             ar(counter).EdgeColor = coWarn;
%             ar(counter).FaceColor = coWarn;
%         end
        
        counter = counter + 1;
        
    end
end


hold off

% set(gca,'visible','off'); % Uncomment to get rid of axis

if nargout == 1
    varargout{1} = f;
elseif nargout == 2
    varargout{1} = f;
    varargout{2} = p;
elseif nargout == 3
    varargout{1} = f;
    varargout{2} = p;
    varargout{3} = ar;
end


end

%% Functions to control ROI behavior on click

function  customWait(hROI)

% Listen for mouse clicks on the ROI
l = addlistener(hROI,'ROIClicked',@clickCallback);

% Block program execution
uiwait;

% Remove listener
delete(l);


end

function clickCallback(~,evt)

% Checks for a double click specifically
if strcmp(evt.SelectionType,'double')
    uiresume;
end

end