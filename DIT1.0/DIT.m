clear all
clc
close all



%%%%%%%%%%%%%%%%%%%  Forwards modeling %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Setting up constants  %%%%%%%%%%%%%%%%%%%

UN=0;
UE=0;
UZ=0;

source=[34.2130,-118.5370]; %Source Coordinates=[lat,lon]
Hyp=[15.00,19.30];	%%% Hypocenter coordinates on the fault=[x,z] (according to SIV standards: http://equake-rc.info/)

strike=122.00; % Fault's Strike
dip=40.00; % Fault's Dip
% rake=101.00; % not usuing a point source
% depth=17.50; % not using a point source
length=18.00; %Fault's Length in km
width=24.00; %Fault's Width in km
n_l=36; %Number of subfaults along the fault's length
n_w=48; %Number of subfaults along the fault's width
Mw=6.71; %Moment Magnitude 
seismic_moment=1.30E+19; %Sesmic moment

open=0.0;

%%%%%%%%%%%%%%%%%%% Reading Source info  %%%%%%%%%%%%%%%%%%%%%

m_true_header=46; %Number of headlines in the file containing subfaults details
m_true_name='s1994NORTHR01ZENG.fsp'; % name of the file containing subfaults details
delimiterIn=' ';
coordname='Coordinates.txt'; % name of the file contaning stations coordinates.

source_details=importdata(m_true_name,delimiterIn,m_true_header);
source_inf=[];
fac=4; %Coarsing Factor in each row (the number of cells are reduced by 4 in each row)

seg_l=(length/n_l)*fac;
seg_w=(width/n_w)*fac;

for jj=(1:(n_w/fac))
	for ii=(1:(n_l/fac))
		comp=[((n_l*fac*(jj-1))+(0*n_l)+(fac*(ii-1))+1),((n_l*fac*(jj-1))+(0*n_l)+(fac*(ii-1))+2),((n_l*fac*(jj-1))+(0*n_l)+(fac*(ii-1))+3),((n_l*fac*(jj-1))+(0*n_l)+(fac*(ii-1))+4),...
            ((n_l*fac*(jj-1))+(1*n_l)+(fac*(ii-1))+1),((n_l*fac*(jj-1))+(1*n_l)+(fac*(ii-1))+2),((n_l*fac*(jj-1))+(1*n_l)+(fac*(ii-1))+3),((n_l*fac*(jj-1))+(1*n_l)+(fac*(ii-1))+4),...
            ((n_l*fac*(jj-1))+(2*n_l)+(fac*(ii-1))+1),((n_l*fac*(jj-1))+(2*n_l)+(fac*(ii-1))+2),((n_l*fac*(jj-1))+(2*n_l)+(fac*(ii-1))+3),((n_l*fac*(jj-1))+(2*n_l)+(fac*(ii-1))+4),...
            ((n_l*fac*(jj-1))+(3*n_l)+(fac*(ii-1))+1),((n_l*fac*(jj-1))+(3*n_l)+(fac*(ii-1))+2),((n_l*fac*(jj-1))+(3*n_l)+(fac*(ii-1))+3),((n_l*fac*(jj-1))+(3*n_l)+(fac*(ii-1))+4)];
		comp_source_slip=(source_details.data(comp(1),6)+source_details.data(comp(2),6)+source_details.data(comp(3),6)+source_details.data(comp(4),6)+source_details.data(comp(5),6)+...
            source_details.data(comp(6),6)+source_details.data(comp(7),6)+source_details.data(comp(8),6)+source_details.data(comp(9),6)+source_details.data(comp(10),6)+...
            source_details.data(comp(11),6)+source_details.data(comp(12),6)+source_details.data(comp(13),6)+source_details.data(comp(14),6)+source_details.data(comp(15),6)+...
            source_details.data(comp(16),6))/16;
		s_lat=source_details.data(comp(1),1);
		s_lon=source_details.data(comp(1),2);
		s_X=source_details.data(comp(1),3);
		s_Y=source_details.data(comp(1),4);
		s_Z=source_details.data(comp(1),5);
		s_slip=comp_source_slip;
		s_rake=source_details.data(comp(1),7);
		source_inf=[source_inf;s_lat,s_lon,s_X,s_Y,abs(s_Z),s_slip,s_rake];
	end
end

n_l=n_l/fac;
n_w=n_w/fac;

%%%%%%%%%%%%%%%%%%% Reading Station Coordinaes  %%%%%%%%%%%%%%

delimiterIn=' ';
headerlineIn=1;
data=importdata(coordname,delimiterIn,headerlineIn);

ind=1;
siz=size(data.data);

def=[];
G=zeros((3*siz(1)),(n_l*n_w));
slip=zeros((n_l*n_w),1);

for sl=1:(n_l*n_w)
    slip(sl,1)=source_inf(sl,6);
end

%%%% Rewrite the program to calculate okada85 with unit slip and store real
%%%% slip in another vector, then multiply them to find real deformation
%%%% (real data). use okada85 output to form G matrix

ok_slip=1.0;

for line=1:siz(1)
	for ns=1:(n_l*n_w)
		st_E=data.data(line,3);
		st_N=data.data(line,2);
		deg=bearing(source_inf(ns,1),source_inf(ns,2),st_N,st_E);
		dist=Haversine(source_inf(ns,1),source_inf(ns,2),st_N,st_E);
		E=dist*sin(deg);
		N=dist*cos(deg);
		[UE,UN,UZ]=okada85(E,N,source_inf(ns,5),strike,dip,seg_l,seg_w,source_inf(ns,7),ok_slip,open);
		G(3*line-2,ns)=UE;
        G(3*line-1,ns)=UN;
        G(3*line,ns)=UZ;
	end
end

syn_true=G*slip;

%%% alternative output: http://www.mathworks.com/help/matlab/import_export/write-to-delimited-data-files.html

std_dev=10;
mu=0;  %mean of the errors added
sigma= 20; % standard deviation
er=normrnd(mu,sigma,(size(syn_true)));
syn_noisy=syn_true+(er*1E-2);

%%%%%%%%%%%%%%%%%%%%%%% Different works on G matrix %%%%%%%%%%%%%%%%%%%%%%

mL2=(G'*G)\(G'*syn_noisy); % L2 regression

tol=1E-3; % tolerance
[mL1,p_value]=IRLS(G,syn_noisy,tol); % L1 regression
disp(['L1 Regression P-Value for IRLS method is:',num2str(p_value)])

cov_mL2_k=(sigma^2).* inv(G'*G);
conf_k= 1.96.*(diag(cov_mL2_k).^0.5); % 95 precentile confidence interval, assuming known error dist

res=syn_true-(G*mL2);
st_dv=(norm(res,2))/sqrt((size(G,1))-(size(G,2))); % standard deviation, assuming unkown error dist
cov_mL2_uk= (st_dv^2).* inv(G'*G);
conf_uk= 1.96.*(diag(cov_mL2_uk).^0.5); % 95 percentile confidence interval of unkown error dist

%%%%%%%%%%%%%%%%%%%% Singular Value Decomposition %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Spectrum %%%%%%%%%%%%%%%%%%%%%%%%%%%

[U,S,V]=svd(G,'econ');
index=1:1:min(size(S,1),size(S,2));
disp(['Please use the Spectrum plot and choose tolerance.'])
figure(1)
plot(index,diag(S),'-bo','linewidth',2.1);
legend('Spectrum of coefficient Matrix');
xlabel('i');
ylabel('Singular Values');
title('\bfSpectrum (Singular Values');
tolr=input('Please insert the tolerance:');

%%%%%%%%%%%%%%%%%%%%%  Condition Number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

condn=cond(G);
disp(['System condition number is: ',num2str(condn)])

%%%%%%%%%%%%%%%%%%%%%  Resolution Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rm=V*(V');
Rd=U*(U');
iden=eye(size(Rm));
svdbias=norm((Rm-iden),2);
disp(['Bias Introduced by SVD Method:',num2str(svdbias)])

[nlow,nhigh]=counter(S,tolr);

condtg=diag(S);
condtg=condtg(1)/condtg(nhigh);
disp(['Condition Number after truncating the coeficient matrix is:',num2str(condtg)]);

mapcnt=2;

figure(mapcnt)
rlslip=reshape(slip,[n_l,n_w]);
rlslip=rlslip';
imagesc(rlslip);
colorbar;
colormap('jet'); % Try graph3d for listing of colormaps
title('\bfReal Slip');
mapcnt=mapcnt+1;

% for kk=1:(size(V,1)-nhigh)
%     figure(mapcnt)
%     temp=reshape(V(:,nhigh+kk),[n_l,n_w]);
%     temp=temp';
%     imagesc(temp);
%     colorbar;
%     colormap('jet'); % Try graph3d for listing of colormaps
%     title(['\bfImage of model null space V0,',num2str(nhigh+kk)]);
%     mapcnt=mapcnt+1;
% end

%%% for more info look at: http://stackoverflow.com/questions/3942892/how-do-i-visualize-a-matrix-with-colors-and-values-displayed

disp('Basis vector spaning data null space:');
disp(num2str(U(:,nhigh+1:nhigh+nlow)));

figure(mapcnt)
imagesc(Rm);
colorbar;
colormap('jet'); % Try graph3d for listing of colormaps
title('\bfModel Resolution Matrix');
mapcnt=mapcnt+1;

figure(mapcnt)
diagRm=diag(Rm);
diagRm=reshape(diagRm,[n_l,n_w]);
diagRm=diagRm';
imagesc(diagRm);
colorbar;
colormap('jet'); % Try graph3d for listing of colormaps
title('\bfDiagonal Elements of Model Resolution Matrix');
mapcnt=mapcnt+1;

figure(mapcnt)
imagesc(Rd);
colorbar;
colormap('jet'); % Try graph3d for listing of colormaps
title('\bfData Resolution Matrix');
mapcnt=mapcnt+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TSVD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[UT,ST,VT]=TSVD(G,tolr);

md=V*(inv(S))*(U')*syn_noisy;
mdt=VT*(inv(ST))*(UT')*syn_noisy;
[nnls_m,resnorm,resid]=lsqnonneg(G,syn_noisy);

%%%%%%%%%%%%%%%%%%%% Checkerboard Test  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_chb=zeros(size(G,2),1);

for ii=1:size(u_chb,1)
    if(mod(ii,2)~=0)
        u_chb(ii,1)=1;
    end
end

chb_d=G*u_chb;
mL2_chb=(G'*G)\(G'*chb_d);
[mL1_chb,p_value_chb]=IRLS(G,chb_d,tol);
disp(['L1 Regression P-Value for IRLS method for Checkerboard Test is:',num2str(p_value_chb)])
mL2_chb=reshape(mL2_chb,[n_l,n_w]);
mL2_chb=mL2_chb';
mL1_chb=reshape(mL1_chb,[n_l,n_w]);
mL1_chb=mL1_chb';
md_chb=V*(inv(S))*(U')*chb_d;
mdt_chb=VT*(inv(ST))*(UT')*chb_d;
nnls_m_chb=lsqnonneg(G,chb_d);

%%%%%%%%%%%%%%%%%%%%%%%% Plotting Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(mapcnt)
tempmL2=reshape(mL2,[n_l,n_w]);
tempmL2=tempmL2';
imagesc(tempmL2);
colorbar;
colormap('jet'); % Try graph3d for listing of colormaps
title('\bfLeast Squared Model');
mapcnt=mapcnt+1;

figure(mapcnt)
tempmL1=reshape(mL1,[n_l,n_w]);
tempmL1=tempmL1';
imagesc(tempmL1);
colorbar;
colormap('jet'); % Try graph3d for listing of colormaps
title('\bfL1 Regression Model (IRLS Method)');
mapcnt=mapcnt+1;

figure(mapcnt)
tempmd=reshape(md,[n_l,n_w]);
tempmd=tempmd';
imagesc(tempmd);
colorbar;
colormap('jet'); % Try graph3d for listing of colormaps
title('\bfGneralized Inverse Solution Model (SVD)');
mapcnt=mapcnt+1;

figure(mapcnt)
tempmdt=reshape(mdt,[n_l,n_w]);
tempmdt=tempmdt';
imagesc(tempmdt);
colorbar;
colormap('jet'); % Try graph3d for listing of colormaps
title('\bfTruncated Singular Value Decomposition Solution Model');
mapcnt=mapcnt+1;

figure(mapcnt)
tempnnls_m=reshape(nnls_m,[n_l,n_w]);
tempnnls_m=tempnnls_m';
imagesc(tempnnls_m);
colorbar;
colormap('jet'); % Try graph3d for listing of colormaps
title('\bfNone Negative Least Square Solution Model');
mapcnt=mapcnt+1;

%%%%%%%%%%%%%%%%% Plotting Checker Board Test Results %%%%%%%%%%%%%%%%%%
%%%%% Checker Board Model
figure(mapcnt)
tempchb=reshape(u_chb,[n_l,n_w]);
tempchb=tempchb';
imagesc(tempchb);
colorbar;
colormap('jet'); % Try graph3d for listing of colormaps
title('\bfChecker Board Model');
mapcnt=mapcnt+1;

figure(mapcnt)
imagesc(mL2_chb);
colorbar;
colormap('jet'); % Try graph3d for listing of colormaps
title('\bfLeast Square Solution Model to Checkerboard Test');
mapcnt=mapcnt+1;

figure(mapcnt)
imagesc(mL1_chb);
colorbar;
colormap('jet'); % Try graph3d for listing of colormaps
title('\bfL1 Regression Solution Model to Checkerboard Test');
mapcnt=mapcnt+1;

figure(mapcnt)
tempmd_chb=reshape(md_chb,[n_l,n_w]);
tempmd_chb=tempmd_chb';
imagesc(tempmd_chb);
colorbar;
colormap('jet'); % Try graph3d for listing of colormaps
title('\bfGeneralized Inverse Solution Model to Checkerboard Test');
mapcnt=mapcnt+1;

figure(mapcnt)
tempmdt_chb=reshape(mdt_chb,[n_l,n_w]);
tempmdt_chb=tempmdt_chb';
imagesc(tempmdt_chb);
colorbar;
colormap('jet'); % Try graph3d for listing of colormaps
title('\bfTruncated Singular Value Decomposition Solution Model to Checkerboard Test');
mapcnt=mapcnt+1;

figure(mapcnt)
tempnnls_m_chb=reshape(nnls_m_chb,[n_l,n_w]);
tempnnls_m_chb=tempnnls_m_chb';
imagesc(tempnnls_m_chb);
colorbar;
colormap('jet'); % Try graph3d for listing of colormaps
title('\bfNone Negative Least Square Solution Model to Checkerboard Test');
mapcnt=mapcnt+1;


disp('End of Program')

