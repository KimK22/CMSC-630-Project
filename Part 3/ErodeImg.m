function NewImage = ErodeImg(Image, D)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[W,L] = size(Image);
NewImage = zeros(W,L);

c=floor(size(D,1)/2);
d=floor(size(D,2)/2);

%Intialize a matrix with size of matrix A
for i=1:W-(2*c)
    for j=1:L-(2*d)     
        newMask = Image(i:i+(2*c),j:j+(2*d));      
        Ediff = min(min(newMask - D));
        if Ediff == 0
            NewImage(i,j) = 1;
        end
    end
end

end

