
%Welcome to MEGARA. This is Version 1.6

clear all
clc

%USER INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

objectname=['hd_10167'];


%If objectname is empty then all stellar images are processed
%to run all MARSDEN EXOCOMET stars use 'exocomet'
%to run all MUSICIAN stars use 'musician'
%to run all MUSICIAN TESS stars use 'tess'
%to run a list of targets define a structure array, e.g. objectname.target1='hd_000001';objectname.target2='hd_000002'; ...

%Months/Runs to reduce
reduction_folders={};
%A colon-separated cell array of the directory names to be reduced, e.g. reduction_folders={'Jan 2018'; 'Feb 2018'};
%If reduction_folders is empty then all folders in the directory where pointer_raw_data_directory.m is located are reduced

%Input stellar values (these can be left blank for MUSICIAN or EXOCOMET targets)
teff=[];
logg=[];
vsini=[];

%USER OPTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Choose 1 for on and 0 for off.

star_census=0;% Do you want to redo census of observations? (default 1)
skip_reduction=0;%Do you want to skip reduction? (default 0)
skip_post_reduction=0;%Do you want to skip processing and post-reduction? (default 0)


blue_data_chop_value=300; % number of pixels to cut off the blue end of the spectrum due to incomplete flats (default=500, 300 for Ca H&K)
apply_barycentric_correction=1;%Do you want to apply the barycentric correction? (default 1)

apply_continuum_fit=1;%apply a continuum fit to normalise the data (default=1)
%continuum fitting uses a previous manual fit in the reduced data directory for the star by default. If this does not exist then one of the options below must be set
manual_continuum_fit=0;%Do you want to manually continuum fit? Use this for precision work or non A-F-type stars. This overwrites previous manual fits.
auto_continuum_fit=1;%Do you want to apply an automatic continuum fit? Reccomended only for A-F-type stars
order_merge=1;%Do you want to merge orders? (default 1) Requires apply_continuum_fit=1

overwrite_fulldata=1;%Do you want to redo processing to order-merged spectra? (default 1)
overwrite_post_reduction=1;%Do you want to re-do cross-correlation and related measurements? (default 1)
cross_correlation=1;%Do you want to cross-correlate? Requires order_merge=1 NB this overwrites any previous cross-correlations for this object
rv_measurement=0; %Do you want to measure radial_velocity and vsini? Requires cross_correlation=1
extend_velocities=0; %extends the velocity axes from +/- 200kms to +/- 400kms(default 0)

supress_figures=0; %Do you want to supress figures during running of MEGARA?(default 0) 

ascii_output=0;%Do you want your output files in ascii format?
fits_output=0;%Do you want your output files in fits format?
FAMIAS_output=0;%Do you want your CCFs output ready to load for FAMIAS? (requires overwrite_post_reduction=1 & cross_correlation=1)

%MEGARA CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[final_data,objectname,star_database]=megara_control_file(objectname,overwrite_fulldata,overwrite_post_reduction,blue_data_chop_value,apply_barycentric_correction,apply_continuum_fit,manual_continuum_fit,auto_continuum_fit,order_merge,cross_correlation,teff,logg,vsini,skip_reduction,reduction_folders,ascii_output,fits_output,star_census,FAMIAS_output,rv_measurement,supress_figures,extend_velocities,skip_post_reduction);



