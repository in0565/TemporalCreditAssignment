function axpt=AxesPoint_1raster_new(xnum,ynum)
%% Axes position 만들기
%%   가로 5, 세로 1의 subplot axis를 만들어준다.
% C= axpt(ix, iy, 1, :)
% R= axpt(ix, iy, 2, :)
% x= axpt(ix, iy, 3, :)   

startpoint = [0.3 0.1];
figwidth = [0.60 0.60];
subwidth = [5]; %% 4 2 4 4 4 subfigure의 크기를 다르게 하기 위하여. 
subheight = [0.3 0.3 0.2]; %%총 합이 1이 되도록
interval_x = 0.02; interval_y = 0.05;
interfig_y = 0.05;

dx=[];  dy =[];
dx = (figwidth(1)-(xnum-1)*interval_x)*subwidth/sum(subwidth); cum_dx=cumsum(dx);
dy = (figwidth(2)-(ynum-1)*interval_y-ynum*interfig_y)*subheight;
axpt = zeros(xnum,ynum,3,4);
for ix = 1:xnum
    for iy = 1:ynum
       if ix == 1
        axpt(ix,iy,1,:) = [startpoint(1) startpoint(2)+(ynum-iy)*(interval_y+interfig_y+2*dy(2))+2*dy(2)+2*interfig_y dx(1) dy(1)];   
        axpt(ix,iy,2,:) = [startpoint(1) startpoint(2)+(ynum-iy)*(interval_y+interfig_y+2*dy(2))+dy(2)+interfig_y dx(1) dy(1)];
        axpt(ix,iy,3,:) = [startpoint(1) startpoint(2)+(ynum-iy)*(interval_y+interfig_y+2*dy(2)) dx(1) dy(2)];
       else 
        axpt(ix,iy,1,:) = [startpoint(1)+(ix-1)*(interval_x)+cum_dx(ix-1) startpoint(2)+(ynum-iy)*(interval_y+interfig_y+2*dy(1))+2*dy(2)+2*interfig_y dx(ix) dy(1)];
        axpt(ix,iy,2,:) = [startpoint(1)+(ix-1)*(interval_x)+cum_dx(ix-1) startpoint(2)+(ynum-iy)*(interval_y+interfig_y+2*dy(1))+dy(2)+interfig_y dx(ix) dy(1)];
        axpt(ix,iy,3,:) = [startpoint(1)+(ix-1)*(interval_x)+cum_dx(ix-1) startpoint(2)+(ynum-iy)*(interval_y+interfig_y+2*dy(2)) dx(ix) dy(2)];   
       end
    end
end