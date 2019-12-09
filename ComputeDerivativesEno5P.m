function dux = ComputeDerivativesEno5P(Datau,Datadeltax,Dataxmesh)
% Compute derivative using ENO method
% input: Datau - solutions of PDE where the rows are indexed by x and the
% columns are indexed by t
%
% output: dux - ux where the rows are indexed by x and the columns are indexed by t

[DataLenx , DataLent] = size(Datau);

%% compute finite differences
utem           = Datau;
u2nd           = utem(2:end,:)-utem(1:end-1,:);
u3rd           = u2nd(2:end,:)-u2nd(1:end-1,:);
u4th           = u3rd(2:end,:)-u3rd(1:end-1,:);
u5th           = u4th(2:end,:)-u4th(1:end-1,:);
u2nd           = abs(u2nd);
u3rd           = abs(u3rd);
u4th           = abs(u4th);
u5th           = abs(u5th);

u10            = u2nd(4:end-4,:);
u01            = u2nd(5:end-3,:);
u210           = u3rd(3:end-4,:);
u101           = u3rd(4:end-3,:);
u012           = u3rd(5:end-2,:);
u3210          = u4th(2:end-4,:);
u2101          = u4th(3:end-3,:);
u1012          = u4th(4:end-2,:);
u0123          = u4th(5:end-1,:);
u43210         = u5th(1:end-4,:);
u32101         = u5th(2:end-3,:);
u21012         = u5th(3:end-2,:);
u10123         = u5th(4:end-1,:);
u01234         = u5th(5:end,:);

e2               = ones(2,1);
Diff4Forward1    = spdiags(fdcoeffF(1,Dataxmesh(1),Dataxmesh(1:5)),[0 1 2 3 4],1,DataLenx);
Diff4Forward2    = spdiags(fdcoeffF(1,Dataxmesh(2),Dataxmesh(1:5)),[0 1 2 3 4],1,DataLenx);
mid = fdcoeffF(1,Dataxmesh(3),Dataxmesh(1:5));
mid(:,3) = [];
Diff4Forward34   = spdiags(e2*mid,[0 1 3 4],2,DataLenx);
Diff4Forward     = [Diff4Forward1 ; Diff4Forward2 ; Diff4Forward34];
mid = fdcoeffF(1,Dataxmesh(end-2),Dataxmesh(end-4:end));
mid(:,3) = [];
Diff4Backward34  = spdiags(e2*mid,[DataLenx-6 DataLenx-5  DataLenx-3 DataLenx-2],2,DataLenx);
Diff4Backward1   = spdiags(fdcoeffF(1,Dataxmesh(end),Dataxmesh(end-4:end)),[DataLenx-5 DataLenx-4 DataLenx-3 DataLenx-2 DataLenx-1],1,DataLenx);
Diff4Backward2   = spdiags(fdcoeffF(1,Dataxmesh(end-1),Dataxmesh(end-4:end)),[DataLenx-5 DataLenx-4 DataLenx-3 DataLenx-2 DataLenx-1],1,DataLenx);
Diff4Backward    = [Diff4Backward34 ; Diff4Backward2 ; Diff4Backward1];

% e2               = ones(2,1);
% Diff4Forward12   = spdiags([-25/12*e2  4*e2  -3*e2  4/3*e2  -1/4*e2],[0 1 2 3 4],2,DataLenx)/Datadeltax;
% Diff4Forward34   = spdiags([e2 -8*e2 8*e2 -e2],[0 1 3 4],2,DataLenx)/(12*Datadeltax);
% Diff4Forward     = [Diff4Forward12 ; Diff4Forward34];
% Diff4Backward34  = spdiags([e2 -8*e2 8*e2 -e2],[DataLenx-6 DataLenx-5  DataLenx-3 DataLenx-2],2,DataLenx)/(12*Datadeltax);
% Diff4Backward12  = spdiags([1/4*e2 -4/3*e2  3*e2  -4*e2  25/12*e2],[DataLenx-6 DataLenx-5 DataLenx-4 DataLenx-3 DataLenx-2],2,DataLenx)/Datadeltax;
% Diff4Backward    = [Diff4Backward34 ; Diff4Backward12];

e = ones(DataLenx,1);
Diff43210      = spdiags([1/4*e  -4/3*e   3*e  -4*e   25/12*e],[-4 -3 -2 -1 0],DataLenx,DataLenx)/Datadeltax;
Diff32101      = spdiags([-1/12*e  1/2*e  -3/2*e  5/6*e  1/4*e],[-3 -2 -1 0 1],DataLenx,DataLenx)/Datadeltax;
Diff21012      = spdiags([1/12*e  -8/12*e  8/12*e  -1/12*e],[-2 -1  1 2],DataLenx,DataLenx)/Datadeltax;
Diff10123      = spdiags([-1/4*e  -5/6*e  3/2*e  -1/2*e  1/12*e],[-1 0 1 2 3],DataLenx,DataLenx)/Datadeltax;
Diff01234      = spdiags([-25/12*e  4*e  -3*e  4/3*e  -1/4*e],[0 1 2 3 4],DataLenx,DataLenx)/Datadeltax;


dux            = zeros(DataLenx,DataLent);      % derivative of u - row x - column t
dux(1:4,:)     = Diff4Forward*utem;             % forward scheme for the first 4 points
dux(end-3:end,:)   = Diff4Backward*utem;             % forward scheme for the first 4 points
ut4 = zeros(DataLenx,1);
uvect            = utem(:,1);

% 43210
A = spdiags(ones(5,1)*fdcoeffF(1,Dataxmesh(5),Dataxmesh(1:5)),[-4 -3 -2 -1 0],5,DataLenx);
ut4(1) = A(1,1)*uvect(1);
ut4(2) = [A(2,1) A(2,2)]*uvect(1:2);
ut4(3) = [A(3,1) A(3,2) A(3,3)]*uvect(1:3);
ut4(4) = [A(4,1) A(4,2) A(4,3) A(4,4)]*uvect(1:4);
for i = 5:DataLenx
    A = spdiags(ones(5,1)*fdcoeffF(1,Dataxmesh(i),Dataxmesh(i-4:i)),[-4 -3 -2 -1 0],5,DataLenx);
    ut4(i)=[A(5,1) A(5,2) A(5,3) A(5,4) A(5,5)]*uvect(i-4:i);
end
ut43210 = ut4;

%32101
ut3 = zeros(DataLenx,1);
A = spdiags(ones(5,1)*fdcoeffF(1,Dataxmesh(4),Dataxmesh(1:5)),[-3 -2 -1 0 1],5,DataLenx);
ut3(1) = [A(1,1) A(1,2)]*uvect(1:2);
ut3(2) = [A(2,1) A(2,2) A(2,3)]*uvect(1:3);
ut3(3) = [A(3,1) A(3,2) A(3,3) A(3,4)]*uvect(1:4);
for i = 4:DataLenx-1
    A = spdiags(ones(5,1)*fdcoeffF(1,Dataxmesh(i),Dataxmesh(i-3:i+1)),[-3 -2 -1 0 1],5,DataLenx);
    ut3(i)=[A(4,1) A(4,2) A(4,3) A(4,4) A(4,5)]*uvect(i-3:i+1);
end
ut3(end) = [A(5,2) A(5,3) A(5,4) A(5,5)]*uvect(end-3:end);
ut32101 = ut3;

% %21012 
ut2 = zeros(DataLenx,1);
mid = fdcoeffF(1,Dataxmesh(3),Dataxmesh(1:5));
mid(:,3) = [];
A  = spdiags(ones(5,1)*mid,[-2 -1 1 2],5,DataLenx);
ut2(1) = [A(1,2) A(1,3)]*[uvect(2) uvect(3)]';
ut2(2) = [A(2,1) A(2,3) A(2,4)]*[uvect(1) uvect(3) uvect(4)]';
for i = 3:DataLenx-2
    mid = fdcoeffF(1,Dataxmesh(3),Dataxmesh(1:5));
    mid(:,3) = [];
    A  = spdiags(ones(5,1)*mid,[-2 -1 1 2],5,DataLenx);
    ut2(i)=[A(3,1) A(3,2) A(3,4) A(3,5)]*[uvect(i-2) uvect(i-1) uvect(i+1) uvect(i+2)]';
end
ut2(end-1) = [A(4,2) A(4,3) A(4,5)]*[uvect(end-3) uvect(end-2) uvect(end)]';
ut2(end) = [A(5,3) A(5,4)]*[uvect(end-2) uvect(end-1)]';
ut21012 = ut2;

%10123
ut1 = zeros(DataLenx,1);
A = spdiags(ones(5,1)*fdcoeffF(1,Dataxmesh(2),Dataxmesh(1:5)),[-1 0 1 2 3],5,DataLenx);
ut1(1) = [A(1,1) A(1,2) A(1,3) A(1,4)]*uvect(1:4);

for i = 2:DataLenx-3
    A = spdiags(ones(5,1)*fdcoeffF(1,Dataxmesh(i),Dataxmesh(i-1:i+3)),[-1 0 1 2 3],5,DataLenx);
    ut1(i)=[A(2,1) A(2,2) A(2,3) A(2,4) A(2,5)]*uvect(i-1:i+3);
end
ut1(end-2) = [A(3,2) A(3,3) A(3,4) A(3,5)]*uvect(end-3:end);
ut1(end-1) = [A(4,3) A(4,4) A(4,5)]*uvect(end-2:end);
ut1(end) = [A(5,4) A(5,5)]*uvect(end-1:end);
ut10123 = ut1;

%01234
ut0 = zeros(DataLenx,1);
A = spdiags(ones(5,1)*fdcoeffF(1,Dataxmesh(2),Dataxmesh(1:5)),[0 1 2 3 4],5,DataLenx);

for i = 1:DataLenx-4
    A = spdiags(ones(5,1)*fdcoeffF(1,Dataxmesh(i),Dataxmesh(i:i+4)),[0 1 2 3 4],5,DataLenx);
    ut0(i)=[A(1,1) A(1,2) A(1,3) A(1,4) A(1,5)]*uvect(i:i+4);
end
ut0(end-3) = [A(2,2) A(2,3) A(2,4) A(2,5)]*uvect(end-3:end);
ut0(end-2) = [A(3,3) A(3,4) A(3,5)]*uvect(end-2:end);
ut0(end-1) = [A(4,4) A(4,5)]*uvect(end-1:end);
ut0(end) = [A(5,5)]*uvect(end);
ut01234 = ut0;

for itert = 1 : DataLent
    uvect            = utem(:,itert);
    EnoDiffu = (ut43210(5:end-4)).*(u10(:,itert)<=u01(:,itert)).*(u210(:,itert)<=u101(:,itert)).*(u3210(:,itert)<=u2101(:,itert)).*(u43210(:,itert)<=u32101(:,itert)) ...
        + (ut32101(5:end-4)).*(u10(:,itert)<=u01(:,itert)).*(u210(:,itert)<=u101(:,itert)).*(u3210(:,itert)<=u2101(:,itert)).*(u43210(:,itert)>u32101(:,itert)) ...
        + (ut32101(5:end-4)).*(u10(:,itert)<=u01(:,itert)).*(u210(:,itert)<=u101(:,itert)).*(u3210(:,itert)>u2101(:,itert)).*(u32101(:,itert)<=u21012(:,itert)) ...
        + (ut21012(5:end-4)).*(u10(:,itert)<=u01(:,itert)).*(u210(:,itert)<=u101(:,itert)).*(u3210(:,itert)>u2101(:,itert)).*(u32101(:,itert)>u21012(:,itert)) ...
        + (ut32101(5:end-4)).*(u10(:,itert)<=u01(:,itert)).*(u210(:,itert)>u101(:,itert)).*(u2101(:,itert)<=u1012(:,itert)).*(u32101(:,itert)<=u21012(:,itert)) ...
        + (ut21012(5:end-4)).*(u10(:,itert)<=u01(:,itert)).*(u210(:,itert)>u101(:,itert)).*(u2101(:,itert)<=u1012(:,itert)).*(u32101(:,itert)>u21012(:,itert)) ...
        + (ut21012(5:end-4)).*(u10(:,itert)<=u01(:,itert)).*(u210(:,itert)>u101(:,itert)).*(u2101(:,itert)>u1012(:,itert)).*(u21012(:,itert)<=u10123(:,itert)) ...
        + (ut10123(5:end-4)).*(u10(:,itert)<=u01(:,itert)).*(u210(:,itert)>u101(:,itert)).*(u2101(:,itert)>u1012(:,itert)).*(u21012(:,itert)>u10123(:,itert)) ...
        + (ut32101(5:end-4)).*(u10(:,itert)>u01(:,itert)).*(u101(:,itert)<=u012(:,itert)).*(u2101(:,itert)<=u1012(:,itert)).*(u32101(:,itert)<=u21012(:,itert)) ...
        + (ut21012(5:end-4)).*(u10(:,itert)>u01(:,itert)).*(u101(:,itert)<=u012(:,itert)).*(u2101(:,itert)<=u1012(:,itert)).*(u32101(:,itert)>u21012(:,itert)) ...
        + (ut21012(5:end-4)).*(u10(:,itert)>u01(:,itert)).*(u101(:,itert)<=u012(:,itert)).*(u2101(:,itert)>u1012(:,itert)).*(u21012(:,itert)<=u10123(:,itert)) ...
        + (ut10123(5:end-4)).*(u10(:,itert)>u01(:,itert)).*(u101(:,itert)<=u012(:,itert)).*(u2101(:,itert)>u1012(:,itert)).*(u21012(:,itert)>u10123(:,itert)) ...
        + (ut21012(5:end-4)).*(u10(:,itert)>u01(:,itert)).*(u101(:,itert)>u012(:,itert)).*(u1012(:,itert)<=u0123(:,itert)).*(u21012(:,itert)<=u10123(:,itert)) ...
        + (ut10123(5:end-4)).*(u10(:,itert)>u01(:,itert)).*(u101(:,itert)>u012(:,itert)).*(u1012(:,itert)<=u0123(:,itert)).*(u21012(:,itert)>u10123(:,itert)) ...
        + (ut10123(5:end-4)).*(u10(:,itert)>u01(:,itert)).*(u101(:,itert)>u012(:,itert)).*(u1012(:,itert)>u0123(:,itert)).*(u10123(:,itert)<=u01234(:,itert)) ...
        + (ut01234(5:end-4)).*(u10(:,itert)>u01(:,itert)).*(u101(:,itert)>u012(:,itert)).*(u1012(:,itert)>u0123(:,itert)).*(u10123(:,itert)>u01234(:,itert));
    dux(5:end-4,itert) = EnoDiffu;
end

end