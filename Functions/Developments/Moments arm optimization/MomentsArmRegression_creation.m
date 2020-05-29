function MomentsArmRegression=MomentsArmRegression_creation()

% All coeffs non explained come from 
% Ramsay, J. W., Hunter, B. V., & Gonzalez, R. V. (2009). 
% Muscle moment arm and normalized moment contributions as reference
% data for musculoskeletal elbow and wrist joint models. 
% Journal of Biomechanics, 42(4), 463–473. https://doi.org/10.1016/j.jbiomech.2008.11.035


k=0;

k=k+1;
%% TricepsBrachii1
MomentsArmRegression(k).name='TricepsBrachii1';
MomentsArmRegression(k).regression(1).equation=1;
MomentsArmRegression(k).regression(1).primaryjoint='Ulna';
MomentsArmRegression(k).regression(1).coeffs=[-24.5454 - 8.8691 9.3509 -1.7518 0]';

k=k+1;
%% TricepsBrachii2
MomentsArmRegression(k).name='TricepsBrachii2';
MomentsArmRegression(k).regression(1).equation=1;
MomentsArmRegression(k).regression(1).primaryjoint='Ulna';
MomentsArmRegression(k).regression(1).coeffs=[-24.5454 - 8.8691 9.3509 -1.7518 0]';

k=k+1;
%% BicepsBrachii1
MomentsArmRegression(k).name='BicepsBrachii1';
MomentsArmRegression(k).regression(1).equation=2;
MomentsArmRegression(k).regression(1).primaryjoint='Radius_J1';
MomentsArmRegression(k).regression(1).secondaryjoint='Radius';
MomentsArmRegression(k).regression(1).coeffs=[ 8.4533 36.6147 2.4777 -19.432 2.0571 13.6502 0 0 -5.6172 0 -2.0854 0 0 0 0 0 ]';
MomentsArmRegression(k).regression(2).equation=2;
MomentsArmRegression(k).regression(2).primaryjoint='Radius';
MomentsArmRegression(k).regression(2).secondaryjoint='Radius_J1';
MomentsArmRegression(k).regression(2).coeffs=[ 1.7271 -4.1504 5.3103 -8.197 0.4405 1.0401 2.6866 -5.5137 0.6448 0 -1.0155 0 2.9534 0 -0.5583 0]';


k=k+1;
%% BicepsBrachii2
MomentsArmRegression(k).name='BicepsBrachii2';
MomentsArmRegression(k).regression(1).equation=2;
MomentsArmRegression(k).regression(1).primaryjoint='Radius_J1';
MomentsArmRegression(k).regression(1).secondaryjoint='Radius';
MomentsArmRegression(k).regression(1).coeffs=[ 8.4533 36.6147 2.4777 -19.432 2.0571 13.6502 0 0 -5.6172 0 -2.0854 0 0 0 0 0 ]';
MomentsArmRegression(k).regression(2).equation=2;
MomentsArmRegression(k).regression(2).primaryjoint='Radius';
MomentsArmRegression(k).regression(2).secondaryjoint='Radius_J1';
MomentsArmRegression(k).regression(2).coeffs=[ 1.7271 -4.1504 5.3103 -8.197 0.4405 1.0401 2.6866 -5.5137 0.6448 0 -1.0155 0 2.9534 0 -0.5583 0]';


k=k+1;
%% Brachialis
MomentsArmRegression(k).name='Brachialis';
MomentsArmRegression(k).regression(1).equation=1;
MomentsArmRegression(k).regression(1).primaryjoint='Radius_J1';
MomentsArmRegression(k).regression(1).coeffs=[ 16.1991 -16.1463 24.5512 -6.3335 0 ]';

k=k+1;
%% Brachioradialis
MomentsArmRegression(k).name='Brachioradialis';
MomentsArmRegression(k).regression(1).equation=2;
MomentsArmRegression(k).regression(1).primaryjoint='Radius_J1';
MomentsArmRegression(k).regression(1).secondaryjoint='Radius';
MomentsArmRegression(k).regression(1).coeffs=[15.2564 -11.8355 2.8129 -5.7781 44.8143 0 2.9032 0 0 -13.4956 0 -0.3940 0 0 0 0]';
MomentsArmRegression(k).regression(2).equation=2;
MomentsArmRegression(k).regression(2).primaryjoint='Radius';
MomentsArmRegression(k).regression(2).secondaryjoint='Radius_J1';
MomentsArmRegression(k).regression(2).coeffs=[ 3.8738 -3.1232 -2.3556 6.0596 7.9944 0.0973 -2.9923 -2.0882 -3.9195 -2.2210 -0.1293 0.4969 0.3683 1.3385 0.9909 -0.3279 ]';

k=k+1;
%% PronatorQuadrus
MomentsArmRegression(k).name='PronatorQuadrus';
MomentsArmRegression(k).regression(1).equation=3;
MomentsArmRegression(k).regression(1).primaryjoint='Radius_J1';
MomentsArmRegression(k).regression(1).secondaryjoint='Radius';
MomentsArmRegression(k).regression(1).coeffs=[11.0405 -1.0079 0.3933 -10.4824 -12.1639 -0.4369 36.9174 3.5232 -10.4223 21.2604 -37.2444 10.2666 -11.0060 14.5974 -3.9919 1.7526 -2.0089 0.5460 ]';
MomentsArmRegression(k).regression(2).equation=2;
MomentsArmRegression(k).regression(2).primaryjoint='Radius';
MomentsArmRegression(k).regression(2).secondaryjoint='Radius_J1';
MomentsArmRegression(k).regression(2).coeffs=[ 5.0238 7.6939 -0.2566 0.9826 -3.3182 0 -0.3034  ]';

k=k+1;
%% SupinatorBrevis'=
MomentsArmRegression(k).name='SupinatorBrevis';
MomentsArmRegression(k).regression(1).equation=1;
MomentsArmRegression(k).regression(1).primaryjoint='Radius';
MomentsArmRegression(k).regression(1).coeffs=[ -13.8661 3.4337 ]';








k=k+1;
%% ExtensorDigitorum
MomentsArmRegression(k).name='ExtensorDigitorum';
MomentsArmRegression(k).regression(1).equation=1;
%MomentsArmRegression(k).regression(1).primaryjoint={"WFE"};
MomentsArmRegression(k).regression(1).primaryjoint='Wrist_J1';
MomentsArmRegression(k).regression(2).equation=1;
%MomentsArmRegression(k).regression(2).primaryjoint={"RUD"};
MomentsArmRegression(k).regression(2).primaryjoint='Hand';
MomentsArmRegression(k).regression(1).coeffs=[-14.1276 1.7325]';
MomentsArmRegression(k).regression(2).coeffs=[2.0459 4.5732]';

% MomentsArmRegression(k).regression(3).primaryjoint='Radius_J1';
% % Data for elbowflexion
% % Taken from Gonzalez, R. V., Buchanan, T. S., & Delp, S. L. (1997). How muscle architecture and moment arms affect wrist flexion-extension moments. Journal of Biomechanics, 30(7), 705–712. https://doi.org/10.1016/S0021-9290(97)00015-8
% % From measures in OpenSim
% 
% x=[15, 40, 52.5 , 62, 71, 79 ,86, 93, 100, 106 , 112, 121, 127]*pi/180;
% y=-0.01:0.001:0.002;
% 
% p=polyfit(x,y,2);
% MomentsArmRegression(k).regression(3).equation=1;
% MomentsArmRegression(k).regression(3).coeffs=flip(p)';





k=k+1;
%% FlexorDigitorumSuperior
MomentsArmRegression(k).name='FlexorDigitorumSuperior';
MomentsArmRegression(k).regression(1).equation=1;
MomentsArmRegression(k).regression(1).primaryjoint='WFE';
MomentsArmRegression(k).regression(2).equation=1;
MomentsArmRegression(k).regression(2).primaryjoint='RUD';
MomentsArmRegression(k).regression(1).coeffs=[10.3467 1.0641 1.0495]';
MomentsArmRegression(k).regression(2).coeffs=[1.6252 6.3604]';

%save('MomentsArmRegression.mat','MomentsArmRegression');