%% Setup

clear;clc;close all;

cd('D:\Electrophysiology\Matdelane\');
addpath(genpath('matnwb'));
addpath(genpath('matdelane'));
addpath(genpath('flipv2'));

generateCore();
addpath('fieldtrip');
ft_defaults;
disp("Setup done.");

load("tfrSet\info.mat");

%%

areax = "V1";
spkpath = "spkSet\";
spkfiles = {dir(spkpath).name};
spkfiles = spkfiles(contains(spkfiles, areax));
Nfiles = length(spkfiles);
spkData = cell(Nfiles, 1);

%%

for ik = 1:Nfiles

    tsetx = load(spkpath + spkfiles{ik});
    spkData{ik} = tsetx.xset;
    fprintf(num2str(ik));

end

%%

%%