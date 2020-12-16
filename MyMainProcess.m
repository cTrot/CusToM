clear all
close all
clc
cheminglobal=pwd;
cd(cheminglobal);

%Chargement des tailles
% excel='Tableau_sujets_Civil.xlsx';
excel='Tableau_sujets_Militaire.xlsx';
tab=importdata(excel);

names = tab.textdata(3:end,2);
taille = tab.data(:,5);
masse = tab.data(:,6);
% taille = tab.data(:,7);
% masse = tab.data(:,8);

doss=dir(cheminglobal);
doss={doss([doss.isdir]').name}';
doss(1:2)=[];

SUJETS=doss;
nb_s=length(SUJETS);

if exist(fullfile(pwd, 'meta.mat'),'file')==2
    load('meta')
else
    %% METADATA TREATEMENT
    meta=struct();
    for ii=1:nb_s
        meta(ii).name=SUJETS{ii};
        [~,ind]=intersect(names,SUJETS{ii});
        meta(ii).taille = taille(ind,1);
        meta(ii).masse = masse(ind,1);
        meta(ii).path = fullfile(cheminglobal,SUJETS{ii});
    end
    %% TRIAL LIST
    for ii=1:nb_s
        cd(meta(ii).path)
        c3dFiles = dir ([meta(ii).path filesep '*.c3d']);
        meta(ii).TrialList={c3dFiles.name}';
        meta(ii).TrialList_tasks=setdiff(meta(ii).TrialList,{'ROM.c3d','Statref.c3d'});
    end
    
    %% GENERATE PARAMETERS
    for ii=1:nb_s
        cd(meta(ii).path)
        
        meta(ii).ModelParameters=struct();
        meta(ii).AnalysisParameters=struct();
        
        [meta(ii).ModelParameters,...
            meta(ii).AnalysisParameters] = Parameters(meta(ii));
        
        for jj=1:length(meta(ii).TrialList_tasks)
            meta(ii).ExternalForces(jj).AnalysisParameters = ...
                AnalysisExternalForces( meta(ii).AnalysisParameters, ...
                meta(ii).TrialList_tasks(jj));
            
            if contains(meta(ii).TrialList_tasks(jj),'Marche')
                
                meta(ii).ExternalForcesPrediction(jj).AnalysisParameters =...
                    AnalysisForcePrediction(meta(ii).AnalysisParameters, ...
                    meta(ii).TrialList_tasks(jj));
                
            end
            
        end
        
    end
    save(fullfile(cheminglobal,'meta'),'meta')
    %% Frames of interest - event
    for ii=1:nb_s
        cd(meta(ii).path)
        for jj=1:length(meta(ii).TrialList_tasks)
            
            t= meta(ii).TrialList_tasks(jj);
            
            if contains(t,'Marche') % A faire pour les autres taches.
                
                FootOnPlates = meta(ii).ExternalForces(jj).AnalysisParameters.ExternalForces.Options;
                [meta(ii).FrameOfInterest(jj).Event,...
                    meta(ii).FrameOfInterest(jj).ID]=...
                    DetectOnlyOneFootPhase(t,FootOnPlates);
                
            elseif contains(t,'Course') || contains(t,'Direction')
                
                FootOnPlates = meta(ii).ExternalForces(jj).AnalysisParameters.ExternalForces.Options;
                [meta(ii).FrameOfInterest(jj).Event,...
                    meta(ii).FrameOfInterest(jj).ID]=...
                    DetectRunPhases(t,FootOnPlates);
                
            elseif contains(t,'Saut')
                
                FootOnPlates = meta(ii).ExternalForces(jj).AnalysisParameters.ExternalForces.Options;
                [meta(ii).FrameOfInterest(jj).Event,...
                    meta(ii).FrameOfInterest(jj).ID]=...
                    DetectJumpPhases(t,FootOnPlates);
                
            end
        end
    end
    save(fullfile(cheminglobal,'meta'),'meta')
end
%% TRAITEMENT EMG
% for ii=1:nb_s % Ajouter pour ne pas tenir compte des EMGs mauvais pour les max !!
%     cd(meta(ii).path)
%     EMG_traitementAll(meta(ii));
% end
%% CALIBRATION GEOMETRIQUE + INVERSE KINEMATICS
for ii=1:nb_s
    cd(meta(ii).path)
    
%     CalibrateModelGenerationNum(meta(ii).ModelParameters,meta(ii).AnalysisParameters);
    InverseKinematics(meta(ii).AnalysisParameters);
    
end

%% EXTERNAL FORCE COMPUTATION + INVERSE DYNAMICS


% Chargement de la fonction symbolique de calcul des forces du sac a dos
load('F4.mat')

%Choix de la methode de determination des forces du sac a dos s'appliquant
%sur l'humain
optionRes = 2;

for ii=1:nb_s
    cd(meta(ii).path)
    for jj=1:length(meta(ii).TrialList_tasks)
        if ~contains(meta(ii).TrialList_tasks(jj),'ROM')
            ExternalForcesComputation(meta(ii).ExternalForces(jj).AnalysisParameters);
            
            %Ajout des forces issues du sac a dos
            if contains(meta(ii).TrialList_tasks(jj),'rec') %pour le port de charge reconstruit
                trialName = meta(ii).TrialList_tasks{jj};
                [fextDos,BiomechanicalModel.OsteoArticularModel,posG] = addBackPackForces(BiomechanicalModel.OsteoArticularModel,trialName,optionRes);
            end
            
            InverseDynamics(meta(ii).ExternalForces(jj).AnalysisParameters);
            
%             %Prediction si marche
%             if contains(meta(ii).TrialList_tasks(jj),'Marche')
%                 ExternalForcesComputation(meta(ii).ExternalForcesPrediction(jj).AnalysisParameters,...
%                     meta(ii).ModelParameters);
%                 
%                 load([meta(ii).TrialList_tasks{jj}(1:end-4) '/InverseDynamicsResults']);
%                 ResID(jj).Ext=InverseDynamicsResults;
%                 
%                 InverseDynamics(meta(ii).ExternalForcesPrediction(jj).AnalysisParameters);
%                 load([meta(ii).TrialList_tasks{jj}(1:end-4) '/InverseDynamicsResults']);
%                 ResID(jj).Pred=InverseDynamicsResults;
%             end
        end
    end
%     save('ResID','ResID')
end

%% MUSCLE FORCE COMPUTATION
 [meta(ii).ModelParameters,...
            meta(ii).AnalysisParameters] = Parameters(meta(ii));
        meta(ii).AnalysisParameters.filename=meta(ii).TrialList_tasks([3,7,12,13,19,20]);
for ii=1:nb_s
    cd(meta(ii).path)
    MuscleForcesComputationNum(meta(ii).AnalysisParameters);
end

%% POSTPROCESS PREDICTION MARCHE
% for ii=1:nb_s
%     cd(meta(ii).path)
%     load('BiomechanicalModel.mat')
%     for jj=1:length(meta(ii).TrialList)
%         if contains(meta(ii).TrialList(jj),'Marche')
%             figure
%             PostProcessExternalForceComparison(BiomechanicalModel,...
%                 meta(ii).ExternalForces(jj).AnalysisParameters,...
%                 meta(ii).TrialList{jj}(1:end-4))
%             for kk=1:6
%                 subplot(2,3,kk)
%                 plotVline(224,'k');
%                 plotVline(400,'k');
%                 xlim([224 400])
%             end
%             figure;
%             q = 14:19; p=numSubplots(length(q));
%             for kk=1:length(q)
%                 subplot(p(1),p(2),kk)
%                 plot(ResID(jj).Ext.JointTorques(q(kk),:)','k-','linewidth',2.5); hold on
%                 plot(ResID(jj).Pred.JointTorques(q(kk),:)','b--','linewidth',2.5)
%                 title(BiomechanicalModel.OsteoArticularModel(q(kk)).name)
%                 plotVline(224,'k');
%                 plotVline(400,'k');
%                 xlim([224 400])
%             end
%         end
%     end
% end

%% PostProcess Kinematics residuals

for ii=1:nb_s
    cd(meta(ii).path)
    for jj=1:length(meta(ii).TrialList_tasks)
        load([meta(ii).TrialList_tasks{jj}(1:end-4) '/InverseKinematicsResults']);
        disp([meta(ii).TrialList_tasks{jj}(1:end-4) ' '...
            num2str(mean(mean(InverseKinematicsResults.ReconstructionError)))])
    end
end

%% PostProcess Dynamic residuals

for ii=1:nb_s
    cd(meta(ii).path)
    for jj=1:length(meta(ii).TrialList_tasks)
        load([meta(ii).TrialList_tasks{jj}(1:end-4) '/InverseDynamicsResults']);
        Window = meta(ii).FrameOfInterest(jj).ID(1):meta(ii).FrameOfInterest(jj).ID(2);
        disp([meta(ii).TrialList_tasks{jj}(1:end-4) ', f6dof '...
            num2str(mean(mean(abs(InverseDynamicsResults.DynamicResiduals.f6dof(:,Window)))))])
        disp([meta(ii).TrialList_tasks{jj}(1:end-4) ', t6dof '...
            num2str(mean(mean(abs(InverseDynamicsResults.DynamicResiduals.t6dof(:,Window)))))])
    end
end

%% PostProcess EMG vs Activation

% for ii=1:nb_s
%     cd(meta(ii).path)
%     
%     load('BiomechanicalModel.mat')
%     MusclesNames={BiomechanicalModel.Muscles.name}';
%     for jj=1:length(meta(ii).TrialList_tasks)
%         
%         trial = meta(ii).TrialList_tasks{jj}(1:end-4);
%         load([trial '/Envelope_EMG.mat'])
%         
%         bool=~[Envelope_EMG.isnan]';
%         Envelope = Envelope_EMG(bool);
%         Env_names={Envelope.name}';
%         
%         n_m=length(Env_names);
%         [~,ind1,ind2]=intersect(MusclesNames,Env_names);
%         
%         EMG=[Envelope(ind2).value]';
%         temg=Envelope(1).time;
%         
%         Window = meta(ii).FrameOfInterest(jj).ID(1):meta(ii).FrameOfInterest(jj).ID(2);
%         
%         load([trial '/MuscleForcesComputationResults'])
%         Act=MuscleForcesComputationResults.MuscleActivations(ind1,Window);
%         
%         load([trial '/ExperimentalData.mat'])
%         tact=ExperimentalData.Time(Window);
%         
%         figure
%         p=numSubplots(n_m);
%         for kk=1:n_m
%             subplot(p(1),p(2),kk)
%             PlotEMGvsActivation(temg,EMG(kk,:)',tact,Act(kk,:)')
%             title(Env_names(ind2(kk)))
%         end
%     end
% end

