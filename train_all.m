
clear; clc; close all;
addpath(genpath(cd))

dbs = ['O','C','F','V'];

for i = 1:length(dbs)
fprintf('%d:%s\n',i,dbs(i));

switch dbs(i)
    case 'O'
        [X,Y] = load_data('O');
        [Wlda_O,~,W] = gen_general(X,Y);
        save W_O W;
    case 'M'
        [X,Y] = load_data('M');
        [Wlda_M,~,W] = gen_general(X,Y);
        save W_M W;
    case 'F'
        [X,Y] = load_data('F');
        [Wlda_F,~,W] = gen_general(X,Y);
        save W_F W;
    case 'C'
        [X,Y] = load_data('C');
        [Wlda_C,~,W] = gen_general(X,Y);
        save W_C W;
    case 'V'
        [X,Y] = load_data('V');
        [Wlda_V,~,W] = gen_general(X,Y);
        save W_V W;        
end

end




