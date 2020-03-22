%% calling palm_quickperms to generate the permuted sample
%create permutation order
function PermSamplePALM(EBfile, PSfile)
% Input:
% file EB_abcd_3076.mat for questionnaire based blocks (NS, SB, TW, TR)
% file EB_abcd_3036.mat for genetic based blocks (NS, MZ, DZ, FS, HS, TR)
% Output:
% crosslag_permorder_3076.mat for the permutation samples using the questionnaire-based blocks
% crosslag_permorder_3036.mat for the permutation samples using the genetic-based blocks

clear,clc
load(EBfile);
B(:,5) = 1:1:size(B,1);
[Pset, VG] = palm_quickperms([], B,5001);
%tabulate(VG)
Pset = Pset(:,2:end);
save(PSfile, 'Pset', 'VG');