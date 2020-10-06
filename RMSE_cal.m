function [d,RMSE] = RMSE_cal(Displacementfield, mask)
%RMSE_CAL Calculates RMSE from a selected region of registered images 
% Input parameters: 
        % Displacementfield: displacement field calculated by non-rigid
        % registration
        % mask: the user selected region from non-rigid registered post-ExM
        % image for RMSE calculation
% Output parameters:    
        % d: measurement length
        % RMSE: calculated RMSE (as a function of d)


[Dx,Dy,~] = size(Displacementfield);
if nargin < 2
    mask =  ones(Dx,Dy);
end

% position coordiante matrix
[n1, m1] = meshgrid(1:1:Dy,1:1:Dx);
ax(:,:,1) = m1;
ax(:,:,2) = n1;
% mask to exclude no feature areas
Dis_0 = Displacementfield.*mask;
ax_0 = ax.*mask;
mask2(:,:,1) = mask;
mask2(:,:,2) = mask;
index_1 = mask2(:,:,1)~=0;
index_2 = mask2(:,:,2)~=0;
Dis_0x = Dis_0(:,:,1);
Dis_0y = Dis_0(:,:,2);
ax_0x = ax_0(:,:,1);
ax_0y = ax_0(:,:,2);


% extract feature into n*2 vector
V(:,1) = Dis_0x(index_1); 
V(:,2) = Dis_0y(index_2);     % V is the matix that contains all distortion vectors of extracted features
ax2(:,1) = ax_0x(index_1);
ax2(:,2) = ax_0y(index_2);     % ax2 is the matrix that contains all coordinates of extracted features
[Dxy,~] = size(ax2);
% make matrices that containcircshift of ax2 and V
ax2_0=repmat(ax2,[Dxy,1]);
V_0=repmat(V,[Dxy,1]);
ax2_m=zeros(Dxy*Dxy,2);
V_m=zeros(Dxy*Dxy,2);
for i=1:Dxy:Dxy*Dxy
    ax2_m(i:i+Dxy-1,:)=circshift(ax2_0(i:i+Dxy-1,:),floor(i/Dxy));
    V_m(i:i+Dxy-1,:)=circshift(V_0(i:i+Dxy-1,:),floor(i/Dxy));
end

% calculate RMSE and d
r_V=V_m-V_0;
r_ax=ax2_m-ax2_0;
d=sqrt(sum(r_ax.^2,2));
RMSE=sqrt(sum(r_V.^2,2));

end

