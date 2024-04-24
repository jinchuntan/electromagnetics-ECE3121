function [Vinitdisp]=setVBCs
% Written by Tan Jin Chun
% Last Modified: 4/9/2023

% Run this code first before running laplacesolv.m
% This script sets up a 2D grid (X, Y) representing the spatial coordinates in your domain.
% Four sides of this square are initialized with voltages Va, Vb, Vc, and Vd. These represent the potential on each electrode.
% Function setVBCs sets the voltages on the four sides of the square.
% Assumption Made: The 4 electrodes together makes up a square

% Note Vinitdisp is flipped so that positions in the array
% correspond to positions in the contour plot grid. 

global X Y Vinitdisp;
global Va Vb Vc Vd;
step = 4/20;
[X,Y] = meshgrid(-2:step:2);

% Assuming that the 4 electrodes are supplied 120 V of AC Voltage
Va = 120;
Vb = 0;
Vc = 120;
Vd = 120;

% Pre-assign array for efficiciency allocation purpose
Vinit=zeros(length(X),length(Y));

% Setting boundary values
Vinit(1,:)=Va;
Vinit(:,length(Y))=Vb;
Vinit(length(X),:)=Vc;
Vinit(:,1)=Vd;
Vinitdisp = flipud(Vinit);