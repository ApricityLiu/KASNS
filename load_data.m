function [X,Y] = load_data(db)

fprintf('dataset is : %s\n',db);

if (db == 'O')
    data = load('Olivetti');
    x = data.fea;
    y = double((data.gnd)');
    if(min(y)==0)
        y=y+1;
    end
elseif (db == 'M')
    data = load('MPEG7');
    x = data.fea;
    y = double(data.gnd);
    if(min(y)==0)
        y=y+1;
    end
elseif (db == 'C')
    data = load('COIL20');
    x = data.fea;
    y = data.gnd;
    if(min(y)==0)
        y=y+1;
    end
elseif (db == 'F')
    data = load('fashion_MNIST');
    x = data.data;
    y = data.real_label;
    if(min(y)==0)
        y=y+1;
    end
    Fcluster = length(unique(y));
    xn = []; yn = [];
    for i=1:Fcluster
        ind = find(y==i);
        xn = [xn;x(ind(1:100),:)];
        yn = [yn;y(ind(1:100),:)];
    end
    clear x y;
    x = xn;
    y = yn;
elseif (db == 'V')
    data = load('vehicle');
    xm = double(data.fea);
    x = zeros(size(xm,1),32*32);
    for iv=1:size(xm,1)
       mid = xm(iv,:);
       x(iv,:) = mid;
    end
    y = zeros(size(xm,1),1);
    ym = double(data.gnd);
    valym = unique(ym);
    for jv=1:length(valym)
        ind = find(ym==valym(jv));
       y(ind,:) = jv; 
    end
    Vcluster = length(unique(y));
    xn = []; yn = [];
    for i=1:Vcluster
       indv = find(y==i);
       xn = [xn;x(indv(1:100),:)];
       yn = [yn;y(indv(1:100),:)];
    end
    clear x y;
    x = xn;
    y = yn;
end

[X,Y] = dataprocess(x,y);

end

function [feanew,gndnew] = dataprocess(fea,gnd)

orid = size(fea,2);
oripix = sqrt(orid);

fea_res = zeros(size(fea,1),40*40);
for i=1:size(fea,1)
   data = fea(i,:);
   data = reshape(data,oripix,oripix);
   data = imresize(data,[40,40]);
   fea_res(i,:) = reshape(data,1,40*40);
end
fea = fea_res;

% H = eye(size(fea,1))-(1/(size(fea,1)))*ones(size(fea,1));
% fea = H*fea;

for i = 1:size(fea,2)
   dim = fea(:,i);
   diff = max(dim)-min(dim);
   diff(diff==0) = inf;
   dim = (dim-min(dim))./diff;
   fea(:,i) = dim;
end

feanew = fea;
gndnew = gnd;

end
