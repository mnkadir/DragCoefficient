%%

FN=[1:501];
CD=[];
for ii=1:length(FN);
    iFN=FN(ii);
    FileName= ['fine' num2str(iFN) '.csv'];
    %% Data load

    % data=load('PP800.99.csv');   %sets variable
    data = csvread(FileName,1,0);

    %% Generating Grid for cylinder
    v=0.0008; %Velocity in the run
    n=100;
    dz=0.1;
    xc=0;
    yc=0;
    r=0.05;
    dth=360/n;
    theta=0:dth:360-dth;
    uA=(2*pi*r/n)*dz;

    [X,Y,Z] = cylinder(r,n);
    xx=X(1,1:n)+xc;
    yy=Y(1,1:n)+yc;
    zz=-0.05:0.1:0.05;

    %% Calculating Drag force

    P_drag=zeros(size(data));



        x=data(:,5);
        y=data(:,6);
        z=data(:,7);
        p=data(:,1);

        index1=find(x>=-0.06 & x<=0.06);
        index2=find(y>=-0.06 & y<=0.06);
        index3=find(z>=-0.05& z<=0.05);

        ind=intersect(index1,index2);
        ind=intersect(ind,index3);

        X1=x(ind);
        Y1=y(ind);
        Z1=z(ind);
        P1=p(ind);
        F1 = scatteredInterpolant(X1,Y1,Z1,P1);

        pp=zeros(length(zz),length(xx));
        fx=zeros(size(pp));
        fy=zeros(size(pp));

        for i=1:length(zz)
            h=zz(i)*ones(size(xx));
            temp=F1(xx,yy,h);
            pp(i,:)=temp;
            fx(i,:)=-uA*temp.*cosd(theta);
            fy(i,:)=-uA*temp.*sind(theta);
        end

        Fx=sum(fx,'All'); % contribute to drag
        Fy=sum(fy,'All');

        P_drag=Fx;


    %% Calculating Drag co-efficient
    ii
    iCd=((P_drag*2)/(v*v*2*r*dz))
    CD(ii)=iCd;

end
