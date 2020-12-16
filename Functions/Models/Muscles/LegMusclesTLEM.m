function [Muscles]=LegMusclesTLEM(Muscles,Signe)
% Definition of the leg muscle model
%   This model contains 166 muscles, actuating the hip, the knee and the
%   ankle joint
%
%   Based on:
%	-	V. Carbone et al., �TLEM 2.0 - A comprehensive musculoskeletal geometry 
%	dataset for subject-specific modeling of lower extremity,� 
%	J. Biomech., vol. 48, no. 5, pp. 734�741, 2015.
%
%   INPUT
%   - Muscles: set of muscles (see the Documentation for the structure)
%   - Signe: side of the leg model ('R' for right side or 'L' for left side)
%   OUTPUT
%   - Muscles: new set of muscles (see the Documentation for the structure)
%________________________________________________________
%
% Licence
% Toolbox distributed under GPL 3.0 Licence
%________________________________________________________
%
% Authors : Antoine Muller, Charles Pontonnier, Pierre Puchaud and
% Georges Dumont
%________________________________________________________

if strcmp(Signe,'Right')
    Signe = 'R';
else
    Signe = 'L';
end

s=cell(0);

s=[s;{[Signe 'AdductorBrevisDistal1'],45.0409,11.101,[],0.016738,0,{['AdductorBrevisDistal1Origin1' Signe  'Pelvis'];['AdductorBrevisDistal1Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorBrevisDistal2'],45.0409,11.101,[],0.022453,0,{['AdductorBrevisDistal2Origin1' Signe  'Pelvis'];['AdductorBrevisDistal2Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorBrevisMid1'],47.2993,10.571,[],0.013458,0,{['AdductorBrevisMid1Origin1' Signe  'Pelvis'];['AdductorBrevisMid1Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorBrevisMid2'],47.2993,10.571,[],0.0178,0,{['AdductorBrevisMid2Origin1' Signe  'Pelvis'];['AdductorBrevisMid2Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorBrevisProximal1'],48.8606,10.2332,[],0.012184,0,{['AdductorBrevisProximal1Origin1' Signe  'Pelvis'];['AdductorBrevisProximal1Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorBrevisProximal2'],48.8606,10.2332,[],0.013922,0,{['AdductorBrevisProximal2Origin1' Signe  'Pelvis'];['AdductorBrevisProximal2Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorLongus1'],55.4615,9.9168,[],0.041484,0,{['AdductorLongus1Origin1' Signe  'Pelvis'];['AdductorLongus1Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorLongus2'],55.4615,9.9168,[],0.053513,0,{['AdductorLongus2Origin1' Signe  'Pelvis'];['AdductorLongus2Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorLongus3'],55.4615,9.9168,[],0.065974,0,{['AdductorLongus3Origin1' Signe  'Pelvis'];['AdductorLongus3Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorLongus4'],55.4615,9.9168,[],0.07842,0,{['AdductorLongus4Origin1' Signe  'Pelvis'];['AdductorLongus4Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorLongus5'],55.4615,9.9168,[],0.090865,0,{['AdductorLongus5Origin1' Signe  'Pelvis'];['AdductorLongus5Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorLongus6'],55.4615,9.9168,[],0.10307,0,{['AdductorLongus6Origin1' Signe  'Pelvis'];['AdductorLongus6Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorMagnusDistal1'],296.3277,10.1802,[],0.21875,0,{['AdductorMagnusDistal1Origin1' Signe  'Pelvis'];['AdductorMagnusDistal1Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorMagnusDistal2'],296.3277,10.1802,[],0.21491,0,{['AdductorMagnusDistal2Origin1' Signe  'Pelvis'];['AdductorMagnusDistal2Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorMagnusDistal3'],296.3277,10.1802,[],0.21072,0,{['AdductorMagnusDistal3Origin1' Signe  'Pelvis'];['AdductorMagnusDistal3Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorMagnusMid1'],89.5221,9.4949,[],0.022218,0,{['AdductorMagnusMid1Origin1' Signe  'Pelvis'];['AdductorMagnusMid1Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorMagnusMid2'],89.5221,9.4949,[],0.040267,0,{['AdductorMagnusMid2Origin1' Signe  'Pelvis'];['AdductorMagnusMid2Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorMagnusMid3'],89.5221,9.4949,[],0.061532,0,{['AdductorMagnusMid3Origin1' Signe  'Pelvis'];['AdductorMagnusMid3Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorMagnusMid4'],89.5221,9.4949,[],0.084092,0,{['AdductorMagnusMid4Origin1' Signe  'Pelvis'];['AdductorMagnusMid4Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorMagnusMid5'],89.5221,9.4949,[],0.10844,0,{['AdductorMagnusMid5Origin1' Signe  'Pelvis'];['AdductorMagnusMid5Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorMagnusMid6'],89.5221,9.4949,[],0.13752,0,{['AdductorMagnusMid6Origin1' Signe  'Pelvis'];['AdductorMagnusMid6Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorMagnusProximal1'],39.2265,9.5599,[],0.043587,0,{['AdductorMagnusProximal1Origin1' Signe  'Pelvis'];['AdductorMagnusProximal1Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorMagnusProximal2'],39.2265,9.5599,[],0.03741,0,{['AdductorMagnusProximal2Origin1' Signe  'Pelvis'];['AdductorMagnusProximal2Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorMagnusProximal3'],39.2265,9.5599,[],0.032548,0,{['AdductorMagnusProximal3Origin1' Signe  'Pelvis'];['AdductorMagnusProximal3Insertion1' Signe  'Thigh']},{};...
[Signe 'AdductorMagnusProximal4'],39.2265,9.5599,[],0.028946,0,{['AdductorMagnusProximal4Origin1' Signe  'Pelvis'];['AdductorMagnusProximal4Insertion1' Signe  'Thigh']},{};...
[Signe 'BicepsFemorisCaputBreve1'],102.1598,9.7886,[],0.14753,0,{['BicepsFemorisCaputBreve1Origin1' Signe  'Thigh'];['BicepsFemorisCaputBreve1Insertion1' Signe  'Shank']},{};...
[Signe 'BicepsFemorisCaputBreve2'],102.1598,9.7886,[],0.096394,0,{['BicepsFemorisCaputBreve2Origin1' Signe  'Thigh'];['BicepsFemorisCaputBreve2Insertion1' Signe  'Shank']},{};...
[Signe 'BicepsFemorisCaputBreve3'],102.1598,9.7886,[],0.037927,0,{['BicepsFemorisCaputBreve3Origin1' Signe  'Thigh'];['BicepsFemorisCaputBreve3Insertion1' Signe  'Shank']},{};...
[Signe 'BicepsFemorisCaputLongum1'],669.2476,8.2929,[],0.33659,0.42993,{['BicepsFemorisCaputLongum1Origin1' Signe  'Pelvis'];['BicepsFemorisCaputLongum1Insertion1' Signe  'Shank']},{};...
[Signe 'ExtensorDigitorumLongus1'],85.8709,5.0949,[],0.42653,0.19692,{['ExtensorDigitorumLongus1Origin1' Signe  'Shank'];['ExtensorDigitorumLongus1Via2' Signe  'Shank'];['ExtensorDigitorumLongus1Via3' Signe  'Shank'];['ExtensorDigitorumLongus1Via4' Signe  'Shank'];['ExtensorDigitorumLongus1Via5' Signe  'Shank'];['ExtensorDigitorumLongus1Via1' Signe  'Foot'];['ExtensorDigitorumLongus1Via2' Signe  'Foot'];['ExtensorDigitorumLongus1Insertion3' Signe  'Foot']},{};...
[Signe 'ExtensorDigitorumLongus2'],85.8709,5.0949,[],0.33666,0.19692,{['ExtensorDigitorumLongus2Origin1' Signe  'Shank'];['ExtensorDigitorumLongus2Via2' Signe  'Shank'];['ExtensorDigitorumLongus2Via3' Signe  'Shank'];['ExtensorDigitorumLongus2Via4' Signe  'Shank'];['ExtensorDigitorumLongus2Via5' Signe  'Shank'];['ExtensorDigitorumLongus2Via1' Signe  'Foot'];['ExtensorDigitorumLongus2Via2' Signe  'Foot'];['ExtensorDigitorumLongus2Insertion3' Signe  'Foot']},{};...
[Signe 'ExtensorDigitorumLongus3'],85.8709,5.0949,[],0.2137,0.19692,{['ExtensorDigitorumLongus3Origin1' Signe  'Shank'];['ExtensorDigitorumLongus3Via2' Signe  'Shank'];['ExtensorDigitorumLongus3Via3' Signe  'Shank'];['ExtensorDigitorumLongus3Via4' Signe  'Shank'];['ExtensorDigitorumLongus3Via5' Signe  'Shank'];['ExtensorDigitorumLongus3Via1' Signe  'Foot'];['ExtensorDigitorumLongus3Via2' Signe  'Foot'];['ExtensorDigitorumLongus3Insertion3' Signe  'Foot']},{};...
[Signe 'ExtensorDigitorumLongus4'],85.8709,5.0949,[],0.20441,0.19692,{['ExtensorDigitorumLongus4Origin1' Signe  'Shank'];['ExtensorDigitorumLongus4Via2' Signe  'Shank'];['ExtensorDigitorumLongus4Via3' Signe  'Shank'];['ExtensorDigitorumLongus4Via4' Signe  'Shank'];['ExtensorDigitorumLongus4Via5' Signe  'Shank'];['ExtensorDigitorumLongus4Via1' Signe  'Foot'];['ExtensorDigitorumLongus4Via2' Signe  'Foot'];['ExtensorDigitorumLongus4Insertion3' Signe  'Foot']},{};...
[Signe 'ExtensorHallucisLongus1'],52.4473,5.0845,[],0.33215,0.30868,{['ExtensorHallucisLongus1Origin1' Signe  'Shank'];['ExtensorHallucisLongus1Via2' Signe  'Shank'];['ExtensorHallucisLongus1Via3' Signe  'Shank'];['ExtensorHallucisLongus1Via4' Signe  'Shank'];['ExtensorHallucisLongus1Via1' Signe  'Foot'];['ExtensorHallucisLongus1Insertion2' Signe  'Foot']},{};...
[Signe 'ExtensorHallucisLongus2'],52.4473,5.0845,[],0.25964,0.30868,{['ExtensorHallucisLongus2Origin1' Signe  'Shank'];['ExtensorHallucisLongus2Via2' Signe  'Shank'];['ExtensorHallucisLongus2Via3' Signe  'Shank'];['ExtensorHallucisLongus2Via4' Signe  'Shank'];['ExtensorHallucisLongus2Via1' Signe  'Foot'];['ExtensorHallucisLongus2Insertion2' Signe  'Foot']},{};...
[Signe 'ExtensorHallucisLongus3'],52.4473,5.0845,[],0.20282,0.30868,{['ExtensorHallucisLongus3Origin1' Signe  'Shank'];['ExtensorHallucisLongus3Via2' Signe  'Shank'];['ExtensorHallucisLongus3Via3' Signe  'Shank'];['ExtensorHallucisLongus3Via4' Signe  'Shank'];['ExtensorHallucisLongus3Via1' Signe  'Foot'];['ExtensorHallucisLongus3Insertion2' Signe  'Foot']},{};...
[Signe 'FlexorDigitorumLongus1'],95.8555,3.2601,[],0.42003,0.3636,{['FlexorDigitorumLongus1Origin1' Signe  'Shank'];['FlexorDigitorumLongus1Via2' Signe  'Shank'];['FlexorDigitorumLongus1Via3' Signe  'Shank'];['FlexorDigitorumLongus1Via4' Signe  'Shank'];['FlexorDigitorumLongus1Via5' Signe  'Shank'];['FlexorDigitorumLongus1Via6' Signe  'Shank'];['FlexorDigitorumLongus1Via7' Signe  'Shank'];['FlexorDigitorumLongus1Via1' Signe  'Foot'];['FlexorDigitorumLongus1Via2' Signe  'Foot'];['FlexorDigitorumLongus1Insertion3' Signe  'Foot']},{};...
[Signe 'FlexorDigitorumLongus2'],95.8555,3.2601,[],0.36944,0.3636,{['FlexorDigitorumLongus2Origin1' Signe  'Shank'];['FlexorDigitorumLongus2Via2' Signe  'Shank'];['FlexorDigitorumLongus2Via3' Signe  'Shank'];['FlexorDigitorumLongus2Via4' Signe  'Shank'];['FlexorDigitorumLongus2Via5' Signe  'Shank'];['FlexorDigitorumLongus2Via6' Signe  'Shank'];['FlexorDigitorumLongus2Via7' Signe  'Shank'];['FlexorDigitorumLongus2Via1' Signe  'Foot'];['FlexorDigitorumLongus2Via2' Signe  'Foot'];['FlexorDigitorumLongus2Insertion3' Signe  'Foot']},{};...
[Signe 'FlexorDigitorumLongus3'],95.8555,3.2601,[],0.31293,0.3636,{['FlexorDigitorumLongus3Origin1' Signe  'Shank'];['FlexorDigitorumLongus3Via2' Signe  'Shank'];['FlexorDigitorumLongus3Via3' Signe  'Shank'];['FlexorDigitorumLongus3Via4' Signe  'Shank'];['FlexorDigitorumLongus3Via5' Signe  'Shank'];['FlexorDigitorumLongus3Via6' Signe  'Shank'];['FlexorDigitorumLongus3Via7' Signe  'Shank'];['FlexorDigitorumLongus3Via1' Signe  'Foot'];['FlexorDigitorumLongus3Via2' Signe  'Foot'];['FlexorDigitorumLongus3Insertion3' Signe  'Foot']},{};...
[Signe 'FlexorDigitorumLongus4'],95.8555,3.2601,[],0.30196,0.3636,{['FlexorDigitorumLongus4Origin1' Signe  'Shank'];['FlexorDigitorumLongus4Via2' Signe  'Shank'];['FlexorDigitorumLongus4Via3' Signe  'Shank'];['FlexorDigitorumLongus4Via4' Signe  'Shank'];['FlexorDigitorumLongus4Via5' Signe  'Shank'];['FlexorDigitorumLongus4Via6' Signe  'Shank'];['FlexorDigitorumLongus4Via7' Signe  'Shank'];['FlexorDigitorumLongus4Via1' Signe  'Foot'];['FlexorDigitorumLongus4Via2' Signe  'Foot'];['FlexorDigitorumLongus4Insertion3' Signe  'Foot']},{};...
[Signe 'FlexorHallucisLongus1'],230.4957,2.1692,[],0.39452,0.42361,{['FlexorHallucisLongus1Origin1' Signe  'Shank'];['FlexorHallucisLongus1Via2' Signe  'Shank'];['FlexorHallucisLongus1Via3' Signe  'Shank'];['FlexorHallucisLongus1Via4' Signe  'Shank'];['FlexorHallucisLongus1Via5' Signe  'Shank'];['FlexorHallucisLongus1Via6' Signe  'Shank'];['FlexorHallucisLongus1Via7' Signe  'Shank'];['FlexorHallucisLongus1Via1' Signe  'Foot'];['FlexorHallucisLongus1Insertion2' Signe  'Foot']},{};...
[Signe 'FlexorHallucisLongus2'],230.4957,2.1692,[],0.3181,0.42361,{['FlexorHallucisLongus2Origin1' Signe  'Shank'];['FlexorHallucisLongus2Via2' Signe  'Shank'];['FlexorHallucisLongus2Via3' Signe  'Shank'];['FlexorHallucisLongus2Via4' Signe  'Shank'];['FlexorHallucisLongus2Via5' Signe  'Shank'];['FlexorHallucisLongus2Via6' Signe  'Shank'];['FlexorHallucisLongus2Via7' Signe  'Shank'];['FlexorHallucisLongus2Via1' Signe  'Foot'];['FlexorHallucisLongus2Insertion2' Signe  'Foot']},{};...
[Signe 'FlexorHallucisLongus3'],230.4957,2.1692,[],0.24821,0.42361,{['FlexorHallucisLongus3Origin1' Signe  'Shank'];['FlexorHallucisLongus3Via2' Signe  'Shank'];['FlexorHallucisLongus3Via3' Signe  'Shank'];['FlexorHallucisLongus3Via4' Signe  'Shank'];['FlexorHallucisLongus3Via5' Signe  'Shank'];['FlexorHallucisLongus3Via6' Signe  'Shank'];['FlexorHallucisLongus3Via7' Signe  'Shank'];['FlexorHallucisLongus3Via1' Signe  'Foot'];['FlexorHallucisLongus3Insertion2' Signe  'Foot']},{};...
[Signe 'GastrocnemiusLateralis1'],546.425,4.9412,[],0.36536,0.36694,{['GastrocnemiusLateralis1Origin1' Signe  'Thigh'];['GastrocnemiusLateralis1Via1' Signe  'Shank'];['GastrocnemiusLateralis1Insertion1' Signe  'Foot']},{};...
[Signe 'GastrocnemiusMedialis1'],1006.9473,5.3131,[],0.36804,0.17914,{['GastrocnemiusMedialis1Origin1' Signe  'Thigh'];['GastrocnemiusMedialis1Via2' Signe  'Thigh'];['GastrocnemiusMedialis1Via1' Signe  'Shank'];['GastrocnemiusMedialis1Via2' Signe  'Shank'];['GastrocnemiusMedialis1Via3' Signe  'Shank'];['GastrocnemiusMedialis1Via4' Signe  'Shank'];['GastrocnemiusMedialis1Via5' Signe  'Shank'];['GastrocnemiusMedialis1Via6' Signe  'Shank'];['GastrocnemiusMedialis1Insertion1' Signe  'Foot']},{};...
[Signe 'GemellusInferior1'],110.762,2.2571,[],0.023292,0,{['GemellusInferior1Origin1' Signe  'Pelvis'];['GemellusInferior1Insertion1' Signe  'Thigh']},{};...
[Signe 'GemellusSuperior1'],110.762,2.2571,[],0.053737,0,{['GemellusSuperior1Origin1' Signe  'Pelvis'];['GemellusSuperior1Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMaximusInferior1'],164.534,16.0048,[],0.056302,0,{['GluteusMaximusInferior1Origin1' Signe  'Pelvis'];['GluteusMaximusInferior1Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMaximusInferior2'],164.534,16.0048,[],0.046707,0,{['GluteusMaximusInferior2Origin1' Signe  'Pelvis'];['GluteusMaximusInferior2Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMaximusInferior3'],164.534,16.0048,[],0.071564,0,{['GluteusMaximusInferior3Origin1' Signe  'Pelvis'];['GluteusMaximusInferior3Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMaximusInferior4'],164.534,16.0048,[],0.058949,0,{['GluteusMaximusInferior4Origin1' Signe  'Pelvis'];['GluteusMaximusInferior4Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMaximusInferior5'],164.534,16.0048,[],0.081813,0,{['GluteusMaximusInferior5Origin1' Signe  'Pelvis'];['GluteusMaximusInferior5Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMaximusInferior6'],164.534,16.0048,[],0.078605,0,{['GluteusMaximusInferior6Origin1' Signe  'Pelvis'];['GluteusMaximusInferior6Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMaximusSuperior1'],85.414,12.6833,[],0.10677,0,{['GluteusMaximusSuperior1Origin1' Signe  'Pelvis'];['GluteusMaximusSuperior1Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMaximusSuperior2'],85.414,12.6833,[],0.10446,0,{['GluteusMaximusSuperior2Origin1' Signe  'Pelvis'];['GluteusMaximusSuperior2Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMaximusSuperior3'],85.414,12.6833,[],0.10123,0,{['GluteusMaximusSuperior3Origin1' Signe  'Pelvis'];['GluteusMaximusSuperior3Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMaximusSuperior4'],85.414,12.6833,[],0.091809,0,{['GluteusMaximusSuperior4Origin1' Signe  'Pelvis'];['GluteusMaximusSuperior4Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMaximusSuperior5'],85.414,12.6833,[],0.092831,0,{['GluteusMaximusSuperior5Origin1' Signe  'Pelvis'];['GluteusMaximusSuperior5Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMaximusSuperior6'],85.414,12.6833,[],0.085671,0,{['GluteusMaximusSuperior6Origin1' Signe  'Pelvis'];['GluteusMaximusSuperior6Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMediusAnterior1'],161.9274,3.8598,[],0.082193,0,{['GluteusMediusAnterior1Origin1' Signe  'Pelvis'];['GluteusMediusAnterior1Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMediusAnterior2'],161.9274,3.8598,[],0.086366,0,{['GluteusMediusAnterior2Origin1' Signe  'Pelvis'];['GluteusMediusAnterior2Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMediusAnterior3'],161.9274,3.8598,[],0.067034,0,{['GluteusMediusAnterior3Origin1' Signe  'Pelvis'];['GluteusMediusAnterior3Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMediusAnterior4'],161.9274,3.8598,[],0.0719,0,{['GluteusMediusAnterior4Origin1' Signe  'Pelvis'];['GluteusMediusAnterior4Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMediusAnterior5'],161.9274,3.8598,[],0.060053,0,{['GluteusMediusAnterior5Origin1' Signe  'Pelvis'];['GluteusMediusAnterior5Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMediusAnterior6'],161.9274,3.8598,[],0.061159,0,{['GluteusMediusAnterior6Origin1' Signe  'Pelvis'];['GluteusMediusAnterior6Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMediusPosterior1'],284.151,4.3991,[],0.11854,0.25816,{['GluteusMediusPosterior1Origin1' Signe  'Pelvis'];['GluteusMediusPosterior1Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMediusPosterior2'],284.151,4.3991,[],0.099324,0.25816,{['GluteusMediusPosterior2Origin1' Signe  'Pelvis'];['GluteusMediusPosterior2Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMediusPosterior3'],284.151,4.3991,[],0.10535,0.25816,{['GluteusMediusPosterior3Origin1' Signe  'Pelvis'];['GluteusMediusPosterior3Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMediusPosterior4'],284.151,4.3991,[],0.081803,0.25816,{['GluteusMediusPosterior4Origin1' Signe  'Pelvis'];['GluteusMediusPosterior4Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMediusPosterior5'],284.151,4.3991,[],0.077415,0.25816,{['GluteusMediusPosterior5Origin1' Signe  'Pelvis'];['GluteusMediusPosterior5Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMediusPosterior6'],284.151,4.3991,[],0.061822,0.25816,{['GluteusMediusPosterior6Origin1' Signe  'Pelvis'];['GluteusMediusPosterior6Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMinimusAnterior1'],230.0578,2.8254,[],0.083413,0,{['GluteusMinimusAnterior1Origin1' Signe  'Pelvis'];['GluteusMinimusAnterior1Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMinimusAnterior2'],230.0578,2.8254,[],0.069072,0,{['GluteusMinimusAnterior2Origin1' Signe  'Pelvis'];['GluteusMinimusAnterior2Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMinimusMid1'],216.832,2.9977,[],0.082376,0,{['GluteusMinimusMid1Origin1' Signe  'Pelvis'];['GluteusMinimusMid1Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMinimusMid2'],216.832,2.9977,[],0.055851,0,{['GluteusMinimusMid2Origin1' Signe  'Pelvis'];['GluteusMinimusMid2Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMinimusPosterior1'],242.149,2.6843,[],0.070035,0,{['GluteusMinimusPosterior1Origin1' Signe  'Pelvis'];['GluteusMinimusPosterior1Insertion1' Signe  'Thigh']},{};...
[Signe 'GluteusMinimusPosterior2'],242.149,2.6843,[],0.047336,0,{['GluteusMinimusPosterior2Origin1' Signe  'Pelvis'];['GluteusMinimusPosterior2Insertion1' Signe  'Thigh']},{};...
[Signe 'Gracilis1'],95.3844,15.2016,[],0.19586,0,{['Gracilis1Origin1' Signe  'Pelvis'];['Gracilis1Via1' Signe  'Shank'];['Gracilis1Via2' Signe  'Shank'];['Gracilis1Via3' Signe  'Shank'];['Gracilis1Insertion4' Signe  'Shank']},{};...
[Signe 'Gracilis2'],95.3844,15.2016,[],0.2515,0,{['Gracilis2Origin1' Signe  'Pelvis'];['Gracilis2Via1' Signe  'Shank'];['Gracilis2Via2' Signe  'Shank'];['Gracilis2Via3' Signe  'Shank'];['Gracilis2Insertion4' Signe  'Shank']},{};...
[Signe 'IliacusLateralis1'],96.8254,7.4877,[],0.078959,0.59941,{['IliacusLateralis1Origin1' Signe  'Pelvis'];['IliacusLateralis1Via2' Signe  'Pelvis'];['IliacusLateralis1Insertion1' Signe  'Thigh']},{};...
[Signe 'IliacusLateralis2'],96.8254,7.4877,[],0.054981,0.59941,{['IliacusLateralis2Origin1' Signe  'Pelvis'];['IliacusLateralis2Via2' Signe  'Pelvis'];['IliacusLateralis2Insertion1' Signe  'Thigh']},{};...
[Signe 'IliacusMedialis1'],86.8708,8.3457,[],0.10659,0,{['IliacusMedialis1Origin1' Signe  'Pelvis'];['IliacusMedialis1Via2' Signe  'Pelvis'];['IliacusMedialis1Insertion1' Signe  'Thigh']},{};...
[Signe 'IliacusMedialis2'],86.8708,8.3457,[],0.098824,0,{['IliacusMedialis2Origin1' Signe  'Pelvis'];['IliacusMedialis2Via2' Signe  'Pelvis'];['IliacusMedialis2Insertion1' Signe  'Thigh']},{};...
[Signe 'IliacusMid1'],171.7073,4.2223,[],0.12797,0,{['IliacusMid1Origin1' Signe  'Pelvis'];['IliacusMid1Via2' Signe  'Pelvis'];['IliacusMid1Insertion1' Signe  'Thigh']},{};...
[Signe 'IliacusMid2'],171.7073,4.2223,[],0.10804,0,{['IliacusMid2Origin1' Signe  'Pelvis'];['IliacusMid2Via2' Signe  'Pelvis'];['IliacusMid2Insertion1' Signe  'Thigh']},{};...
[Signe 'ObturatorExternusInferior1'],63.564,5.8996,[],0.011732,0,{['ObturatorExternusInferior1Origin1' Signe  'Pelvis'];['ObturatorExternusInferior1Via1' Signe  'Thigh'];['ObturatorExternusInferior1Insertion2' Signe  'Thigh']},{};...
[Signe 'ObturatorExternusInferior2'],63.564,5.8996,[],0.020972,0,{['ObturatorExternusInferior2Origin1' Signe  'Pelvis'];['ObturatorExternusInferior2Via1' Signe  'Thigh'];['ObturatorExternusInferior2Insertion2' Signe  'Thigh']},{};...
[Signe 'ObturatorExternusSuperior1'],182.272,2.3774,[],0.063652,0,{['ObturatorExternusSuperior1Origin1' Signe  'Pelvis'];['ObturatorExternusSuperior1Via1' Signe  'Thigh'];['ObturatorExternusSuperior1Insertion2' Signe  'Thigh']},{};...
[Signe 'ObturatorExternusSuperior2'],182.272,2.3774,[],0.072533,0,{['ObturatorExternusSuperior2Origin1' Signe  'Pelvis'];['ObturatorExternusSuperior2Via1' Signe  'Thigh'];['ObturatorExternusSuperior2Insertion2' Signe  'Thigh']},{};...
[Signe 'ObturatorExternusSuperior3'],182.272,2.3774,[],0.082093,0,{['ObturatorExternusSuperior3Origin1' Signe  'Pelvis'];['ObturatorExternusSuperior3Via1' Signe  'Thigh'];['ObturatorExternusSuperior3Insertion2' Signe  'Thigh']},{};...
[Signe 'ObturatorInternus1'],145.532,1.7751,[],0.12122,0,{['ObturatorInternus1Origin1' Signe  'Pelvis'];['ObturatorInternus1Via2' Signe  'Pelvis'];['ObturatorInternus1Insertion1' Signe  'Thigh']},{};...
[Signe 'ObturatorInternus2'],145.532,1.7751,[],0.10557,0,{['ObturatorInternus2Origin1' Signe  'Pelvis'];['ObturatorInternus2Via2' Signe  'Pelvis'];['ObturatorInternus2Insertion1' Signe  'Thigh']},{};...
[Signe 'ObturatorInternus3'],145.532,1.7751,[],0.12557,0,{['ObturatorInternus3Origin1' Signe  'Pelvis'];['ObturatorInternus3Via2' Signe  'Pelvis'];['ObturatorInternus3Insertion1' Signe  'Thigh']},{};...
[Signe 'ObturatorInternus4'],145.532,1.7751,[],0.080203,0,{['ObturatorInternus4Origin1' Signe  'Pelvis'];['ObturatorInternus4Via2' Signe  'Pelvis'];['ObturatorInternus4Insertion1' Signe  'Thigh']},{};...
[Signe 'ObturatorInternus5'],145.532,1.7751,[],0.11925,0,{['ObturatorInternus5Origin1' Signe  'Pelvis'];['ObturatorInternus5Via2' Signe  'Pelvis'];['ObturatorInternus5Insertion1' Signe  'Thigh']},{};...
[Signe 'ObturatorInternus6'],145.532,1.7751,[],0.085628,0,{['ObturatorInternus6Origin1' Signe  'Pelvis'];['ObturatorInternus6Via2' Signe  'Pelvis'];['ObturatorInternus6Insertion1' Signe  'Thigh']},{};...
[Signe 'Pectineus1'],49.4079,9.3608,[],0.003881,0,{['Pectineus1Origin1' Signe  'Pelvis'];['Pectineus1Insertion1' Signe  'Thigh']},{};...
[Signe 'Pectineus2'],49.4079,9.3608,[],0.0085835,0,{['Pectineus2Origin1' Signe  'Pelvis'];['Pectineus2Insertion1' Signe  'Thigh']},{};...
[Signe 'Pectineus3'],49.4079,9.3608,[],0.012277,0,{['Pectineus3Origin1' Signe  'Pelvis'];['Pectineus3Insertion1' Signe  'Thigh']},{};...
[Signe 'Pectineus4'],49.4079,9.3608,[],0.016838,0,{['Pectineus4Origin1' Signe  'Pelvis'];['Pectineus4Insertion1' Signe  'Thigh']},{};...
[Signe 'PeroneusBrevis1'],142.4393,2.2232,[],0.26155,0.446,{['PeroneusBrevis1Origin1' Signe  'Shank'];['PeroneusBrevis1Via2' Signe  'Shank'];['PeroneusBrevis1Via3' Signe  'Shank'];['PeroneusBrevis1Via4' Signe  'Shank'];['PeroneusBrevis1Via1' Signe  'Foot'];['PeroneusBrevis1Via2' Signe  'Foot'];['PeroneusBrevis1Insertion3' Signe  'Foot']},{};...
[Signe 'PeroneusBrevis2'],142.4393,2.2232,[],0.2505,0.446,{['PeroneusBrevis2Origin1' Signe  'Shank'];['PeroneusBrevis2Via2' Signe  'Shank'];['PeroneusBrevis2Via3' Signe  'Shank'];['PeroneusBrevis2Via4' Signe  'Shank'];['PeroneusBrevis2Via1' Signe  'Foot'];['PeroneusBrevis2Via2' Signe  'Foot'];['PeroneusBrevis2Insertion3' Signe  'Foot']},{};...
[Signe 'PeroneusBrevis3'],142.4393,2.2232,[],0.16523,0.446,{['PeroneusBrevis3Origin1' Signe  'Shank'];['PeroneusBrevis3Via2' Signe  'Shank'];['PeroneusBrevis3Via3' Signe  'Shank'];['PeroneusBrevis3Via4' Signe  'Shank'];['PeroneusBrevis3Via1' Signe  'Foot'];['PeroneusBrevis3Via2' Signe  'Foot'];['PeroneusBrevis3Insertion3' Signe  'Foot']},{};...
[Signe 'PeroneusLongus1'],236.3285,2.962,[],0.39684,0.294,{['PeroneusLongus1Origin1' Signe  'Shank'];['PeroneusLongus1Via2' Signe  'Shank'];['PeroneusLongus1Via3' Signe  'Shank'];['PeroneusLongus1Via4' Signe  'Shank'];['PeroneusLongus1Via5' Signe  'Shank'];['PeroneusLongus1Via1' Signe  'Foot'];['PeroneusLongus1Via2' Signe  'Foot'];['PeroneusLongus1Via3' Signe  'Foot'];['PeroneusLongus1Insertion4' Signe  'Foot']},{};...
[Signe 'PeroneusLongus2'],236.3285,2.962,[],0.34032,0.294,{['PeroneusLongus2Origin1' Signe  'Shank'];['PeroneusLongus2Via2' Signe  'Shank'];['PeroneusLongus2Via3' Signe  'Shank'];['PeroneusLongus2Via4' Signe  'Shank'];['PeroneusLongus2Via5' Signe  'Shank'];['PeroneusLongus2Via1' Signe  'Foot'];['PeroneusLongus2Via2' Signe  'Foot'];['PeroneusLongus2Via3' Signe  'Foot'];['PeroneusLongus2Insertion4' Signe  'Foot']},{};...
[Signe 'PeroneusLongus3'],236.3285,2.962,[],0.26875,0.294,{['PeroneusLongus3Origin1' Signe  'Shank'];['PeroneusLongus3Via2' Signe  'Shank'];['PeroneusLongus3Via3' Signe  'Shank'];['PeroneusLongus3Via4' Signe  'Shank'];['PeroneusLongus3Via5' Signe  'Shank'];['PeroneusLongus3Via1' Signe  'Foot'];['PeroneusLongus3Via2' Signe  'Foot'];['PeroneusLongus3Via3' Signe  'Foot'];['PeroneusLongus3Insertion4' Signe  'Foot']},{};...
[Signe 'Piriformis1'],310.6871,4.0233,[],0.087382,0,{['Piriformis1Origin1' Signe  'Pelvis'];['Piriformis1Insertion1' Signe  'Thigh']},{};...
[Signe 'Plantaris1'],58.6392,4.2634,[],0.35838,0,{['Plantaris1Origin1' Signe  'Thigh'];['Plantaris1Insertion1' Signe  'Foot']},{};...
[Signe 'Popliteus1'],147.0918,2.0395,[],0.043274,0,{['Popliteus1Origin1' Signe  'Thigh'];['Popliteus1Via1' Signe  'Shank'];['Popliteus1Insertion2' Signe  'Shank']},{};...
[Signe 'Popliteus2'],147.0918,2.0395,[],0.063865,0,{['Popliteus2Origin1' Signe  'Thigh'];['Popliteus2Via1' Signe  'Shank'];['Popliteus2Insertion2' Signe  'Shank']},{};...
[Signe 'Popliteus3'],147.0918,2.0395,[],0.097843,0,{['Popliteus3Origin1' Signe  'Thigh'];['Popliteus3Via1' Signe  'Shank'];['Popliteus3Insertion2' Signe  'Shank']},{};...
[Signe 'PsoasMajor1'],109.4432,9.1372,[],0.22357,0.28083,{['PsoasMajorT12I_TMVia5' Signe  'Pelvis'];['PsoasMajor1Via1' Signe  'Pelvis'];['PsoasMajor1Insertion1' Signe  'Thigh']},{};...
[Signe 'PsoasMajor2'],109.4432,9.1372,[],0.19218,0.28083,{['PsoasMajor2Origin1' Signe  'Spine*'];['PsoasMajor2Via1' Signe  'Pelvis'];['PsoasMajor2Insertion1' Signe  'Thigh']},{};...
[Signe 'PsoasMajor3'],109.4432,9.1372,[],0.1648,0.28083,{['PsoasMajor3Origin1' Signe  'Spine*'];['PsoasMajor3Via1' Signe  'Pelvis'];['PsoasMajor3Insertion1' Signe  'Thigh']},{};...
[Signe 'PsoasMajor4'],109.4432,9.1372,[],0.13728,0.28083,{['PsoasMajor4Origin1' Signe  'Spine*'];['PsoasMajor4Via1' Signe  'Pelvis'];['PsoasMajor4Insertion1' Signe  'Thigh']},{};...
[Signe 'PsoasMajor5'],109.4432,9.1372,[],0.11063,0.28083,{['PsoasMajor5Origin1' Signe  'Spine*'];['PsoasMajor5Via1' Signe  'Pelvis'];['PsoasMajor5Insertion1' Signe  'Thigh']},{};...
[Signe 'QuadratusFemoris1'],142.763,2.8894,[],0.019905,0,{['QuadratusFemoris1Origin1' Signe  'Pelvis'];['QuadratusFemoris1Insertion1' Signe  'Thigh']},{};...
[Signe 'QuadratusFemoris2'],142.763,2.8894,[],0.024002,0,{['QuadratusFemoris2Origin1' Signe  'Pelvis'];['QuadratusFemoris2Insertion1' Signe  'Thigh']},{};...
[Signe 'QuadratusFemoris3'],142.763,2.8894,[],0.032185,0,{['QuadratusFemoris3Origin1' Signe  'Pelvis'];['QuadratusFemoris3Insertion1' Signe  'Thigh']},{};...
[Signe 'QuadratusFemoris4'],142.763,2.8894,[],0.042174,0,{['QuadratusFemoris4Origin1' Signe  'Pelvis'];['QuadratusFemoris4Insertion1' Signe  'Thigh']},{};...
[Signe 'RectusFemoris1'],390.0439,7.3069,[],0.30049,0.32581,{['RectusFemoris1Origin1' Signe  'Pelvis'];['RectusFemoris1Insertion1' Signe  'Patella']},{};...
[Signe 'RectusFemoris2'],390.0439,7.3069,[],0.29326,0.32581,{['RectusFemoris2Origin1' Signe  'Pelvis'];['RectusFemoris2Insertion1' Signe  'Patella']},{};...
[Signe 'Sartorius1'],149.3788,32.8025,[],0.11608,0,{['Sartorius1Origin1' Signe  'Pelvis'];['Sartorius1Via1' Signe  'Thigh'];['Sartorius1Via2' Signe  'Thigh'];['Sartorius1Via3' Signe  'Thigh'];['Sartorius1Via4' Signe  'Thigh'];['Sartorius1Via5' Signe  'Thigh'];['Sartorius1Via1' Signe  'Shank'];['Sartorius1Via2' Signe  'Shank'];['Sartorius1Via3' Signe  'Shank'];['Sartorius1Insertion4' Signe  'Shank']},{};...
[Signe 'Semimembranosus1'],247.0407,7.826,[],0.33808,0.37805,{['Semimembranosus1Origin1' Signe  'Pelvis'];['Semimembranosus1Insertion1' Signe  'Shank']},{};...
[Signe 'Semimembranosus2'],247.0407,7.826,[],0.32647,0.37805,{['Semimembranosus2Origin1' Signe  'Pelvis'];['Semimembranosus2Insertion1' Signe  'Shank']},{};...
[Signe 'Semimembranosus3'],247.0407,7.826,[],0.31198,0.37805,{['Semimembranosus3Origin1' Signe  'Pelvis'];['Semimembranosus3Insertion1' Signe  'Shank']},{};...
[Signe 'Semitendinosus1'],410.3055,12.9172,[],0.29889,0,{['Semitendinosus1Origin1' Signe  'Pelvis'];['Semitendinosus1Via1' Signe  'Shank'];['Semitendinosus1Insertion2' Signe  'Shank']},{};...
[Signe 'SoleusLateralis1'],1116.3228,2.1798,[],0.31549,0.68767,{['SoleusLateralis1Origin1' Signe  'Shank'];['SoleusLateralis1Insertion1' Signe  'Foot']},{};...
[Signe 'SoleusLateralis2'],1116.3228,2.1798,[],0.26552,0.68767,{['SoleusLateralis2Origin1' Signe  'Shank'];['SoleusLateralis2Insertion1' Signe  'Foot']},{};...
[Signe 'SoleusLateralis3'],1116.3228,2.1798,[],0.19876,0.68767,{['SoleusLateralis3Origin1' Signe  'Shank'];['SoleusLateralis3Insertion1' Signe  'Foot']},{};...
[Signe 'SoleusMedialis1'],653.9587,2.0389,[],0.26296,0.72828,{['SoleusMedialis1Origin1' Signe  'Shank'];['SoleusMedialis1Insertion1' Signe  'Foot']},{};...
[Signe 'SoleusMedialis2'],653.9587,2.0389,[],0.23934,0.72828,{['SoleusMedialis2Origin1' Signe  'Shank'];['SoleusMedialis2Insertion1' Signe  'Foot']},{};...
[Signe 'SoleusMedialis3'],653.9587,2.0389,[],0.21815,0.72828,{['SoleusMedialis3Origin1' Signe  'Shank'];['SoleusMedialis3Insertion1' Signe  'Foot']},{};...
[Signe 'TensorFasciaeLatae1'],98.6829,8.6135,[],0.3599,0,{['TensorFasciaeLatae1Origin1' Signe  'Pelvis'];['TensorFasciaeLatae1Insertion1' Signe  'Shank']},{};...
[Signe 'TensorFasciaeLatae2'],98.6829,8.6135,[],0.35264,0,{['TensorFasciaeLatae2Origin1' Signe  'Pelvis'];['TensorFasciaeLatae2Insertion1' Signe  'Shank']},{};...
[Signe 'TibialisAnterior1'],319.8905,3.9076,[],0.30953,0.20977,{['TibialisAnterior1Origin1' Signe  'Shank'];['TibialisAnterior1Via2' Signe  'Shank'];['TibialisAnterior1Via1' Signe  'Foot'];['TibialisAnterior1Insertion2' Signe  'Foot']},{};...
[Signe 'TibialisAnterior2'],319.8905,3.9076,[],0.25352,0.20977,{['TibialisAnterior2Origin1' Signe  'Shank'];['TibialisAnterior2Via2' Signe  'Shank'];['TibialisAnterior2Via1' Signe  'Foot'];['TibialisAnterior2Insertion2' Signe  'Foot']},{};...
[Signe 'TibialisAnterior3'],319.8905,3.9076,[],0.1719,0.20977,{['TibialisAnterior3Origin1' Signe  'Shank'];['TibialisAnterior3Via2' Signe  'Shank'];['TibialisAnterior3Via1' Signe  'Foot'];['TibialisAnterior3Insertion2' Signe  'Foot']},{};...
[Signe 'TibialisPosteriorLateralis1'],333.0303,2.152,[],0.3238,0.56776,{['TibialisPosteriorLateralis1Origin1' Signe  'Shank'];['TibialisPosteriorLateralis1Via2' Signe  'Shank'];['TibialisPosteriorLateralis1Via3' Signe  'Shank'];['TibialisPosteriorLateralis1Via4' Signe  'Shank'];['TibialisPosteriorLateralis1Via1' Signe  'Foot'];['TibialisPosteriorLateralis1Via2' Signe  'Foot'];['TibialisPosteriorLateralis1Via3' Signe  'Foot'];['TibialisPosteriorLateralis1Insertion4' Signe  'Foot']},{};...
[Signe 'TibialisPosteriorLateralis2'],333.0303,2.152,[],0.25273,0.56776,{['TibialisPosteriorLateralis2Origin1' Signe  'Shank'];['TibialisPosteriorLateralis2Via2' Signe  'Shank'];['TibialisPosteriorLateralis2Via3' Signe  'Shank'];['TibialisPosteriorLateralis2Via4' Signe  'Shank'];['TibialisPosteriorLateralis2Via1' Signe  'Foot'];['TibialisPosteriorLateralis2Via2' Signe  'Foot'];['TibialisPosteriorLateralis2Via3' Signe  'Foot'];['TibialisPosteriorLateralis2Insertion4' Signe  'Foot']},{};...
[Signe 'TibialisPosteriorLateralis3'],333.0303,2.152,[],0.17027,0.56776,{['TibialisPosteriorLateralis3Origin1' Signe  'Shank'];['TibialisPosteriorLateralis3Via2' Signe  'Shank'];['TibialisPosteriorLateralis3Via3' Signe  'Shank'];['TibialisPosteriorLateralis3Via4' Signe  'Shank'];['TibialisPosteriorLateralis3Via1' Signe  'Foot'];['TibialisPosteriorLateralis3Via2' Signe  'Foot'];['TibialisPosteriorLateralis3Via3' Signe  'Foot'];['TibialisPosteriorLateralis3Insertion4' Signe  'Foot']},{};...
[Signe 'TibialisPosteriorMedialis1'],333.0303,2.152,[],0.33055,0.34168,{['TibialisPosteriorMedialis1Origin1' Signe  'Shank'];['TibialisPosteriorMedialis1Via2' Signe  'Shank'];['TibialisPosteriorMedialis1Via3' Signe  'Shank'];['TibialisPosteriorMedialis1Via4' Signe  'Shank'];['TibialisPosteriorMedialis1Via1' Signe  'Foot'];['TibialisPosteriorMedialis1Via2' Signe  'Foot'];['TibialisPosteriorMedialis1Via3' Signe  'Foot'];['TibialisPosteriorMedialis1Insertion4' Signe  'Foot']},{};...
[Signe 'TibialisPosteriorMedialis2'],333.0303,2.152,[],0.25857,0.34168,{['TibialisPosteriorMedialis2Origin1' Signe  'Shank'];['TibialisPosteriorMedialis2Via2' Signe  'Shank'];['TibialisPosteriorMedialis2Via3' Signe  'Shank'];['TibialisPosteriorMedialis2Via4' Signe  'Shank'];['TibialisPosteriorMedialis2Via1' Signe  'Foot'];['TibialisPosteriorMedialis2Via2' Signe  'Foot'];['TibialisPosteriorMedialis2Via3' Signe  'Foot'];['TibialisPosteriorMedialis2Insertion4' Signe  'Foot']},{};...
[Signe 'TibialisPosteriorMedialis3'],333.0303,2.152,[],0.16677,0.34168,{['TibialisPosteriorMedialis3Origin1' Signe  'Shank'];['TibialisPosteriorMedialis3Via2' Signe  'Shank'];['TibialisPosteriorMedialis3Via3' Signe  'Shank'];['TibialisPosteriorMedialis3Via4' Signe  'Shank'];['TibialisPosteriorMedialis3Via1' Signe  'Foot'];['TibialisPosteriorMedialis3Via2' Signe  'Foot'];['TibialisPosteriorMedialis3Via3' Signe  'Foot']},{};...
[Signe 'VastusIntermedius1'],121.3047,6.9384,[],0.1906,0.16486,{['VastusIntermedius1Origin1' Signe  'Thigh'];['VastusIntermedius1Insertion1' Signe  'Patella']},{};...
[Signe 'VastusIntermedius2'],121.3047,6.9384,[],0.18963,0.16486,{['VastusIntermedius2Origin1' Signe  'Thigh'];['VastusIntermedius2Insertion1' Signe  'Patella']},{};...
[Signe 'VastusIntermedius3'],121.3047,6.9384,[],0.1155,0.16486,{['VastusIntermedius3Origin1' Signe  'Thigh'];['VastusIntermedius3Insertion1' Signe  'Patella']},{};...
[Signe 'VastusIntermedius4'],121.3047,6.9384,[],0.11586,0.16486,{['VastusIntermedius4Origin1' Signe  'Thigh'];['VastusIntermedius4Insertion1' Signe  'Patella']},{};...
[Signe 'VastusIntermedius5'],121.3047,6.9384,[],0.03249,0.16486,{['VastusIntermedius5Origin1' Signe  'Thigh'];['VastusIntermedius5Insertion1' Signe  'Patella']},{};...
[Signe 'VastusIntermedius6'],121.3047,6.9384,[],0.033265,0.16486,{['VastusIntermedius6Origin1' Signe  'Thigh'];['VastusIntermedius6Insertion1' Signe  'Patella']},{};...
[Signe 'VastusLateralisInferior1'],211.9371,3.3029,[],0.23118,0,{['VastusLateralisInferior1Origin1' Signe  'Thigh'];['VastusLateralisInferior1Insertion1' Signe  'Patella']},{};...
[Signe 'VastusLateralisInferior2'],211.9371,3.3029,[],0.18008,0,{['VastusLateralisInferior2Origin1' Signe  'Thigh'];['VastusLateralisInferior2Insertion1' Signe  'Patella']},{};...
[Signe 'VastusLateralisInferior3'],211.9371,3.3029,[],0.13869,0,{['VastusLateralisInferior3Origin1' Signe  'Thigh'];['VastusLateralisInferior3Insertion1' Signe  'Patella']},{};...
[Signe 'VastusLateralisInferior4'],211.9371,3.3029,[],0.095187,0,{['VastusLateralisInferior4Origin1' Signe  'Thigh'];['VastusLateralisInferior4Insertion1' Signe  'Patella']},{};...
[Signe 'VastusLateralisInferior5'],211.9371,3.3029,[],0.064617,0,{['VastusLateralisInferior5Origin1' Signe  'Thigh'];['VastusLateralisInferior5Insertion1' Signe  'Patella']},{};...
[Signe 'VastusLateralisInferior6'],211.9371,3.3029,[],0.032468,0,{['VastusLateralisInferior6Origin1' Signe  'Thigh'];['VastusLateralisInferior6Insertion1' Signe  'Patella']},{};...
[Signe 'VastusLateralisSuperior1'],1061.0355,7.7754,[],0.24421,0,{['VastusLateralisSuperior1Origin1' Signe  'Thigh'];['VastusLateralisSuperior1Insertion1' Signe  'Patella']},{};...
[Signe 'VastusLateralisSuperior2'],1061.0355,7.7754,[],0.22506,0,{['VastusLateralisSuperior2Origin1' Signe  'Thigh'];['VastusLateralisSuperior2Insertion1' Signe  'Patella']},{};...
[Signe 'VastusMedialisInferior1'],134.5015,8.5501,[],0.056105,0,{['VastusMedialisInferior1Origin1' Signe  'Thigh'];['VastusMedialisInferior1Insertion1' Signe  'Patella']},{};...
[Signe 'VastusMedialisInferior2'],134.5015,8.5501,[],0.026143,0,{['VastusMedialisInferior2Origin1' Signe  'Thigh'];['VastusMedialisInferior2Insertion1' Signe  'Patella']},{};...
[Signe 'VastusMedialisMid1'],276.1659,7.9662,[],0.11057,0,{['VastusMedialisMid1Origin1' Signe  'Thigh'];['VastusMedialisMid1Insertion1' Signe  'Patella']},{};...
[Signe 'VastusMedialisMid2'],276.1659,7.9662,[],0.080481,0,{['VastusMedialisMid2Origin1' Signe  'Thigh'];['VastusMedialisMid2Insertion1' Signe  'Patella']},{};...
[Signe 'VastusMedialisSuperior1'],143.5881,9.7501,[],0.20255,0,{['VastusMedialisSuperior1Origin1' Signe  'Thigh'];['VastusMedialisSuperior1Insertion1' Signe  'Patella']},{};...
[Signe 'VastusMedialisSuperior2'],143.5881,9.7501,[],0.19456,0,{['VastusMedialisSuperior2Origin1' Signe  'Thigh'];['VastusMedialisSuperior2Insertion1' Signe  'Patella']},{};...
[Signe 'VastusMedialisSuperior3'],143.5881,9.7501,[],0.18049,0,{['VastusMedialisSuperior3Origin1' Signe  'Thigh'];['VastusMedialisSuperior3Insertion1' Signe  'Patella']},{};...
[Signe 'VastusMedialisSuperior4'],143.5881,9.7501,[],0.15846,0,{['VastusMedialisSuperior4Origin1' Signe  'Thigh'];['VastusMedialisSuperior4Insertion1' Signe  'Patella']},{};...
    }];


% Cr�ation de la structure
Muscles=[Muscles;struct('name',{s{:,1}}','f0',{s{:,2}}','l0',{s{:,3}}',...
    'Kt',{s{:,4}}','ls',{s{:,5}}','alpha0',{s{:,6}}','path',{s{:,7}}','wrap',{s{:,8}}')]; %#ok<CCAT1>

end