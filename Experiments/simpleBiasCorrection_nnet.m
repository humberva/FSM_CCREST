function [y1] = simpleBiasCorrection_nnet(x1)
%MYNEURALNETWORKFUNCTION neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 10-May-2017 17:54:22.
%
% [y1] = myNeuralNetworkFunction(x1) takes these arguments:
%   x = 4xQ matrix, input #1
% and returns:
%   y = 1xQ matrix, output #1
% where Q is the number of samples.
%
% X = [log(morphdata(qbias(qc_idx,1),2)), log(climdata(qbias(qc_idx,1),7)), climdata(qbias(qc_idx,1),8), sqrt(surfdata(qbias(qc_idx,1),6))];
%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1_xoffset = [0.0353671438372913;4.80138642963982;-1.54142;0];
x1_step1_gain = [0.134645615919686;0.551137917750831;0.0816225918274564;2.40317100951024];
x1_step1_ymin = -1;

% Layer 1
b1 = [-2.763483693281550746;2.1867187739925251222;2.4732755104101373789;1.123092601114394018;0.077316361930475779873;-0.23647942533859289682;0.58572106065454532242;0.36415495870192426953;1.6214361316828180826;1.5217909434809264724];
IW1_1 = [1.1132433312709679729 -1.9780838557092805807 2.2803626006295898421 -0.37932328003817655793;-0.29526290071003757909 -3.2016794363573026772 2.733920133457097279 -1.527993465618578961;-1.6606734355089769473 -2.6528743047210725514 -1.8025603690423328551 1.0470140051879175402;-1.6887459899976664524 1.659929346090086133 0.74444546796691724033 0.15925924199601793063;0.77075560999391246053 -0.65929829862555267006 0.32094689086515493237 0.15454964988021377414;-1.2758794854450672407 1.8811020093121246788 -1.2040472017245895664 0.23716594800273715915;-1.1427335937119500464 -1.4093505445012566302 -0.085179620964410882045 2.8851019725653994641;0.85550778710376529368 0.89912252862034591772 0.14205350678966532918 -0.77478782064305218302;0.4133575927796483418 -1.6418555697722987397 -0.7751427278218763206 -0.11520569512292938574;1.7708071232769042602 -0.079338959668013755988 -1.3938195133960924466 -0.85695822346580385886];

% Layer 2
b2 = 0.57410838240160355639;
LW2_1 = [0.60498123468445841322 0.27801371889843751184 0.19211851713988978285 0.0011472381040324511337 -1.1662600436121373093 -0.30187600687147841949 -0.065879406313503904857 -0.63591239313989467163 -0.50639623766408004002 0.17283257199911314372];

% Output 1
y1_step1_ymin = -1;
y1_step1_gain = 0.111040173952678;
y1_step1_xoffset = -15.1284471345125;

% ===== SIMULATION ========

% Dimensions
Q = size(x1,2); % samples

% Input 1
xp1 = mapminmax_apply(x1,x1_step1_gain,x1_step1_xoffset,x1_step1_ymin);

% Layer 1
a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*xp1);

% Layer 2
a2 = repmat(b2,1,Q) + LW2_1*a1;

% Output 1
y1 = mapminmax_reverse(a2,y1_step1_gain,y1_step1_xoffset,y1_step1_ymin);
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings_gain,settings_xoffset,settings_ymin)
y = bsxfun(@minus,x,settings_xoffset);
y = bsxfun(@times,y,settings_gain);
y = bsxfun(@plus,y,settings_ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n)
a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings_gain,settings_xoffset,settings_ymin)
x = bsxfun(@minus,y,settings_ymin);
x = bsxfun(@rdivide,x,settings_gain);
x = bsxfun(@plus,x,settings_xoffset);
end