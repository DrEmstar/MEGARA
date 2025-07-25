
%Manual CCF
clear all
load final_data_hd_150549.mat
fullwave=finaldata.fullwave;
fullint=finaldata.fullint;
jd=finaldata.jd;

%Put star values here
objectname=['hd_150549'];
teff=10244;
logg=3.53;
vsini=65.7;
sigm=0.5;wstart=3800;wstop=8000;

weak_limit1=1; %lines < XX mA in strength removed from ccf (default 1)
weak_limit2=5; %lines < YY mA in strength removed from delta function ccf (default 5)

%Cross-correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[synthwave, synthint, wavelengths, ele, elem, widths] = make_synth_spectrum(teff,logg,vsini,sigm,wstart,wstop);
w=wavelengths;
wi=widths;
%cut out Hydrogen and telluric regions from delta function list. The below are the regions to keep.

% %general list for template CCF
 S1=w>4260 & w<4320;
 S2=w>4380 & w<4830;
 S3=w>4895 & w<5865;
 S4=w>5982 & w<6274;
 S5=w>6326 & w<6465;
 S6=w>6600 & w<6858;
 S7=w>7415 & w<7500;

%example of a different region list to cross-corelate
% % S1=w>3985 & w<4030;
% S2=w>4038 & w<4090;
% S3=w>4117 & w<4324;
% S4=w>4356 & w<4798;
% S5=w>4807 & w<4845;
% S6=w>4877 & w<5342;
% S7=w>5345 & w<5387;
% S8=w>5390 & w<5450;
% S9=w>5454 & w<5467;
% S10=w>5472 & w<5515;
% S11=w>5518 & w<5514;
% S12=w>5620 & w<6270;
% S13=w>6325 & w<6522;
% S14=w>6600 & w<6865;
% S15=w>6970 & w<7160;
% S16=w>7375 & w<7590;
% S17=w>7715 & w<8100;
% S18=w>8380 & w<8700;

W=[];
Wi=[];
Ele=[];
list_regions=who('S*');
for n=1:numel(list_regions)
    W=cat(1,W,w(eval(list_regions{n})));
    Wi=cat(1,Wi,wi(eval(list_regions{n})));
    Ele=cat(1,Ele,elem(eval(list_regions{n})));
end

%remove very weak lines
TSS=Wi<weak_limit1;
W(TSS)=[];
Wi(TSS)=[];
Ele(TSS)=[];

%to select only lines of a certain species e.g. FeI, comment out this section if not needed
specie='SiII';
W2=[];
Wi2=[];
for n=1:numel(Ele)
    test_specie=Ele{n};
    if strcmpi(specie,test_specie)
        W2=cat(1,W2,W(n));
        Wi2=cat(1,Wi2,Wi(n));
    end
end
W=W2;
Wi=Wi2;

%do a ccf with the standard synthetic spectrum 

%Find Nans and set to 1
G=isnan(fullint);
fullint(G)=1;

[v,vel,LSD,CCF,dfn_auto,wcw] = do_ccf_from_spec_wavelengths_and_widths(fullwave,fullint,W,Wi,400);

%getting weights from a noise estimate and producing a summed spectrum
intfilt=medfilt1(fullint(:,1:10000)',11)';
resid=fullint(:,1:10000)-intfilt;
weight=1./std(resid,[],2);

G=vel>-200 & vel<200;
U=vel(G);
LSDy=LSD(:,G);

%plot the LSDs
figure
plot(U,LSDy)
eval(['title(''CCF ' objectname(4:end) ''')'])

%coarse_ccf.velocity=U;
%coarse_ccf.intensity=LSDy;
%coarse_ccf.jd=jd;
%eval(['save(''ccf_' objectname '.mat'',''coarse_ccf'')'])
mkdir("Manual_CCF_FAMIAS_files_" + specie)
cd("Manual_CCF_FAMIAS_files_" + specie)
eval(['write_files_for_FAMIAS(jd,U,LSDy,weight,''times_' objectname(4:end) '.txt'',''' objectname(4:end) '_'')'])

