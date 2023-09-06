fname='3';% image label sequence
data=importdata([fname,'-SR-568.txt']);% load supper resolution localization information of channel 1
X1=data.data(:,5);
Y1=data.data(:,6);
Z1=data.data(:,7);
XYZ=unique([X1,Y1,Z1],'rows');
XYZ(isnan(XYZ(:,2)),:)=[];
X1=XYZ(:,1);Y1=XYZ(:,2);Z1=XYZ(:,3);
[v c]=voronoin([X1,Y1,Z1]);
Vol1=zeros(length(c),1);
for i=1:length(c)
    try
        [K av]=convhull([v(c{i},1),v(c{i},2),v(c{i},3)]);
        Vol1(i)=av;
    end
end
[K cv]=convhull(X1,Y1,Z1);
mVol1=cv/length(X1);
Vol1(Vol1==0)=max(Vol1);
IN=Vol1(Vol1~=0)<mVol1;
DT1=delaunayTriangulation([X1,Y1,Z1]);
clear FOid1
FOid1=cell(length(X1),1);
tic
for i=1:length(DT1.ConnectivityList)
    FOid1{DT1.ConnectivityList(i,1)}=[FOid1{DT1.ConnectivityList(i,1)},DT1.ConnectivityList(i,[2,3,4])];
    FOid1{DT1.ConnectivityList(i,2)}=[FOid1{DT1.ConnectivityList(i,2)},DT1.ConnectivityList(i,[1,3,4])];
    FOid1{DT1.ConnectivityList(i,3)}=[FOid1{DT1.ConnectivityList(i,3)},DT1.ConnectivityList(i,[1,2,4])];
    FOid1{DT1.ConnectivityList(i,4)}=[FOid1{DT1.ConnectivityList(i,4)},DT1.ConnectivityList(i,[1,2,3])];
    if mod(i,10000)==0
        i
    end
end
toc
for i=1:length(X1)
    FOid1{i}=unique(FOid1{i});
end

data=importdata([fname,'-SR-647.txt']);% load supper resolution localization information of channel 2
X2=data.data(:,5);
Y2=data.data(:,6);
Z2=data.data(:,7);
XYZ=unique([X2,Y2,Z2],'rows');
X2=XYZ(:,1);Y2=XYZ(:,2);Z2=XYZ(:,3);
[v c]=voronoin([X2,Y2,Z2]);
Vol2=zeros(length(c),1);
for i=1:length(c)
    try
        [K av]=convhull([v(c{i},1),v(c{i},2),v(c{i},3)]);
        Vol2(i)=av;
    end
end
[K cv]=convhull(X2,Y2,Z2);
mVol2=cv/length(X2);
Vol2(Vol2==0)=max(Vol2);
IN=Vol2(Vol2~=0)<mVol2;
DT2=delaunayTriangulation([X2,Y2,Z2]);
clear FOid2
FOid2=cell(length(X2),1);
tic
for i=1:length(DT2.ConnectivityList)
    FOid2{DT2.ConnectivityList(i,1)}=[FOid2{DT2.ConnectivityList(i,1)},DT2.ConnectivityList(i,[2,3,4])];
    FOid2{DT2.ConnectivityList(i,2)}=[FOid2{DT2.ConnectivityList(i,2)},DT2.ConnectivityList(i,[1,3,4])];
    FOid2{DT2.ConnectivityList(i,3)}=[FOid2{DT2.ConnectivityList(i,3)},DT2.ConnectivityList(i,[1,2,4])];
    FOid2{DT2.ConnectivityList(i,4)}=[FOid2{DT2.ConnectivityList(i,4)},DT2.ConnectivityList(i,[1,2,3])];
    if mod(i,10000)==0
        i
    end
end
toc
for i=1:length(X2)
    FOid2{i}=unique(FOid2{i});
end
Vol11=zeros(length(X1),1);
Vol21=zeros(length(X2),1);
for i=1:length(X1)
    Vol11(i)=mean(Vol1([i,FOid1{i}]));
end
for i=1:length(X2)
    Vol21(i)=mean(Vol2([i,FOid2{i}]));
end
save([fname,'-FOid.mat'],'FOid1','DT1','mVol1','Vol1','FOid2','DT2','mVol2','Vol2','Vol11','Vol21','X1','Y1','Z1','X2','Y2','Z2')