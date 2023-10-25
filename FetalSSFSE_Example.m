clear all;
% In this example we simulate Philip's T2-weighted single-shot fast spin
% echo MRI scans of a simulated moving fetal brain. In the future we may
% expand it to multiple vendors. 
% The example below simulates a single-shot fast spin echo sequence on a 33
% weeks old fetal brain. To use this simulator, realistic RF pulses need to
% be used with their appropriate amplitude and duration. 
%%
% Obtained RF pulses from a vendor MRI simulator for a 2.5mm slice.
load("Pulses.mat")
exc.rf = excPulse;
exc.amp = 5.855e-3;
exc.T = 2.893e-3;
exc.rewinder=true;
exc.G = 7.51;
exc.delz = 7.5e-3;
refoc.rf = refPulse;
refoc.amp = [13.671e-3,12.152e-3];
refoc.T = 1.504e-3;
refoc.rewinder=false;
refoc.G = 6.26;
refoc.delz = 7.5e-3;
nAngles = 62;
%%
delz = exc.delz*1000;
NbSlices = 90;
SliceThickness = 2.5; % mm
SliceGap = -1.25; % mm
SliceWindow = delz;
ImRes = 1.25; % mm iso
Acq = 288; % square matrix
%%
TSE_Factor = 99; % Make sure the TSE factor...
FOV = [350,350,NbSlices*ImRes]; % ... and FOV (and resolution) correspond to each other.
B0 = 1.5; 
TE = 180; % Echo time
ESP = 6.1; % Echo spacing
ACF = 2; % Acceleration factor
std_noise = 100; % We tried different stds to get this number as we scale the intensities.
orientation = 1; % Sagittal, 2 = Coronal, 3 = Axial
GA = 33; % Weeks
PhantomRes = 0.8*ones(1,3); %1mm iso resolution of the atlas
MotionLevel = 2; % Moderate Motion:
% Motion levels: 0:4 motion levels: No, Light, Moderate, Severe and Extreme motion.
%%
Fetal_Brain_model_path = 'ExamplePhantom/';
%%
GetCorruptedStacks(Fetal_Brain_model_path,GA,exc,refoc,nAngles,NbSlices, ...
                    SliceWindow,ImRes,Acq,TSE_Factor,FOV,B0,TE,ESP,ACF,MotionLevel,PhantomRes, ...
                    SliceThickness,SliceGap,std_noise,orientation,"GA33_Test");
