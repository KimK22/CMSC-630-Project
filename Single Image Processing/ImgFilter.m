function [NewImage, Hnew] = ImgFilter(Image,color, weight, scale, type)
%ImgFilter: Process Image through a filter given weight inputs from user
ImageCopy = Image(:,:,color);
[W, L] = size(ImageCopy); % width and length of image
inputs = []; % empty vector
Psum = 0; % sum of pixels
maxInt = 255;
NewImage = ImageCopy; % copy original Image

% Linear Filter
if type == 1
    for u = 2:W-1
        for v = 2:L-1
            Psum = 0; % sum of average
            for i = -1:1
                for j = -1:1
                    Psum = Psum + ImageCopy(u+i, v+j)*weight(i+2, j+2)/scale;
                end
            end
            NewImage(u,v) = round(Psum);
        end
    end


% Median Filter
elseif type == 2
    for u = 2:W-1
        for v = 2:L-1
            inputs = [];
            for i = -1:1
                for j = -1:1
                    for val = 1:weight(i+2, j+2)   
                        inputs = [inputs, ImageCopy(u+i,v+j)];
                    end
                end
            end
        NewImage(u,v) = median(inputs);       
        end
    end
end

% Histogram Calculations
for i = 0:maxInt
    Hnew(i+1) = sum(sum(NewImage == i));
    pdf(i+1) = Hnew(i+1) / (W*L); 
        % probability density (# pixels of intensity H / total # pixels)
end

end
