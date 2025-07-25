function extracted = extract_order_no_background(points,ofit,width,img,widthfactor)

%get the size of the image
szimg=size(img);
xax=points;

%define y-position for order
ypositions=round(ofit);

sxax=zeros(size(points));
syax=zeros(numel(points),numel( ceil(-width*widthfactor):floor(width*widthfactor) ));
sdata=zeros(size(syax));
sypositions=sxax;
t=sxax;
cnt=0;

for xpos=1:numel(xax)
    if ceil(ypositions(xpos)-width*widthfactor) > 1 && floor(ypositions(xpos)+width*widthfactor) < szimg(1)
        t(xpos)=1;
        cnt = cnt+1;
        sxax(cnt) = xax(xpos);
        syax(cnt,:) = ceil(ypositions(xpos)-width*widthfactor):floor(ypositions(xpos)+width*widthfactor);
        sdata(cnt,:) = img(ceil(ypositions(xpos)-width*widthfactor):floor(ypositions(xpos)+width*widthfactor), xax(xpos));
        sypositions(cnt)=ofit(xpos);
    end
end

t=round(t)==0;
sxax(t)=[];
sypositions(t)=[];
syax(t,:)=[];
sdata(t,:)=[];

extracted.xax=sxax;
extracted.yax=syax;
extracted.data=sdata;
extracted.ypositions=sypositions;
