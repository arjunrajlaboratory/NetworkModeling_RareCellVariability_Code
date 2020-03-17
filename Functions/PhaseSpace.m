load('/Volumes/MELANOMA/Data/Data1000')
I = [26 92 133 183 544 702 915 968];

xlin = linspace(0.001,0.1,10);
ylin = linspace(0.1,1,10);

countx = 1;
M = zeros(10);

for ixlin = 1:length(xlin)-1
    countx = countx + 1;
    county = 1;
    for iylin = 1:length(ylin)-1
        county = county + 1;
        for idata = I
            if xlin(ixlin)<=Data1000(idata,3) && Data1000(idata,3)<xlin(ixlin+1) 
                if ylin(iylin)<=Data1000(idata,5) && Data1000(idata,5)<ylin(iylin+1)
                    M(county,countx) = M(county,countx) +1;
                end
            end
        end
    end
end

M = flipud(M);
% Vq = interp2(M,10);


%// Define integer grid of coordinates for the above data
[X,Y] = meshgrid(1:size(M,2), 1:size(M,1));

%// Define a finer grid of points
[X2,Y2] = meshgrid(1:0.01:size(M,2), 1:0.01:size(M,1));

%// Interpolate the data and show the output
outData = interp2(X, Y, M, X2, Y2, 'linear');
imagesc(outData);
c = gray;
c = flipud(c);
colormap(c)
colorbar


% imagesc(M)
% hold on
% plot(Data1000(:,3),Data1000(:,5),'.','MarkerSize',10,'Color',[51, 51, 51]./255)
% hold on
% plot(Data1000(I,3),Data1000(I,5),'.','MarkerSize',20,'Color',[255,127,42]./255)
% c = gray;
% c = flipud(c);
% colormap(c)

% figure
% [X,Y] = meshgrid(xlin,ylin);
% contourf(X,Y,M,10,'Linecolor','none')
% % hold on
% % plot(Data1000(:,3),Data1000(:,5),'.','MarkerSize',10,'Color',[51, 51, 51]./255)
% % hold on
% % plot(Data1000(I,3),Data1000(I,5),'.','MarkerSize',20,'Color',[255,127,42]./255)
% c = gray;
% c = flipud(c);
% colormap(c)
% xlim([0,0.1])
% ylim([0,1])

%%

load('Data1000')
I = [26 92 133 183 544 702 915 968];

xlin = linspace(0.001,0.1,10);
ylin = linspace(0.01,0.1,10);

countx = 1;
M = zeros(10);

for ixlin = 1:length(xlin)-1
    countx = countx + 1;
    county = 1;
    for iylin = 1:length(ylin)-1
        county = county + 1;
        for idata = I
            if xlin(ixlin)<=Data1000(idata,3) && Data1000(idata,3)<xlin(ixlin+1) 
                if ylin(iylin)<=Data1000(idata,6) && Data1000(idata,6)<ylin(iylin+1)
                    M(county,countx) = M(county,countx) +1;
                end
            end
        end
    end
end

M = flipud(M);
% Vq = interp2(M,10);


%// Define integer grid of coordinates for the above data
[X,Y] = meshgrid(1:size(M,2), 1:size(M,1));

%// Define a finer grid of points
[X2,Y2] = meshgrid(1:0.01:size(M,2), 1:0.01:size(M,1));

%// Interpolate the data and show the output
outData = interp2(X, Y, M, X2, Y2, 'linear');
imagesc(outData);
c = gray;
c = flipud(c);
colormap(c)
colorbar

% M = flipud(M);
% Vq = interp2(M,10);

% map = [100, 100, 100;
%     255,127,42]./255;
% 
% imagesc(Vq)
% 
% colormap('gray');

% figure
% [X,Y] = meshgrid(xlin,ylin);
% contourf(X,Y,M,10,'Linecolor','none')
% hold on 
% plot(Data1000(:,3),Data1000(:,6),'.','MarkerSize',10,'Color',[51, 51, 51]./255)
% hold on
% plot(Data1000(I,3),Data1000(I,6),'.','MarkerSize',20,'Color',[255,127,42]./255)
% c = gray;
% c = flipud(c);
% colormap(c)
%%

load('Data1000')
I = [26 92 133 183 544 702 915 968];

xlin = linspace(0.01,0.1,10);
ylin = linspace(0.1,1,10);

countx = 1;
M = zeros(10);

for ixlin = 1:length(xlin)-1
    countx = countx + 1;
    county = 1;
    for iylin = 1:length(ylin)-1
        county = county + 1;
        for idata = I
            if xlin(ixlin)<=Data1000(idata,6) && Data1000(idata,6)<xlin(ixlin+1) 
                if ylin(iylin)<=Data1000(idata,5) && Data1000(idata,5)<ylin(iylin+1)
                    M(county,countx) = M(county,countx) +1;
                end
            end
        end
    end
end

M = flipud(M);
% Vq = interp2(M,10);


%// Define integer grid of coordinates for the above data
[X,Y] = meshgrid(1:size(M,2), 1:size(M,1));

%// Define a finer grid of points
[X2,Y2] = meshgrid(1:0.01:size(M,2), 1:0.01:size(M,1));

%// Interpolate the data and show the output
outData = interp2(X, Y, M, X2, Y2, 'linear');
imagesc(outData);
c = gray;
c = flipud(c);
colormap(c)
colorbar

% M = flipud(M);
% Vq = interp2(M,10);

% map = [100, 100, 100;
%     255,127,42]./255;
% 
% imagesc(Vq)
% 
% colormap('gray');

% figure
% [X,Y] = meshgrid(xlin,ylin);
% contourf(X,Y,M,10,'Linecolor','none')
% c = gray;
% c = flipud(c);
% colormap(c)