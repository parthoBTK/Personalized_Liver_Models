function [CondiModel,expressionRxns_GIMME,expressionRxns_EFLUX,genes_rxns_asso_GIMME,genes_rxns_asso_Eflux] = ConditionSpecificEfluxGIMME(refModel,genes,expression_all,grpIND,nom,rxn_cutoff_GIMME)

disp('Starting Model contextualization ...')
disp('GIMME Threshold / Cutoff:')
rxn_cutoff_GIMME

%%% This function sequentially constrian the reference model to GIMME and
%%% E-flux

% G = refModel.genes;
% E = zeros(length(refModel.genes), size(grpIND,2));
% 
% for uu = 1:length(G)
%     uu
%     idpp = find(ismember(genes,G(uu)));
%     if(idpp>0)
%      E(uu,:) = expression_all(idpp(1),grpIND);
%     end
%     
% end

%E(find(ismember(genes,refModel.genes)),grpIND) = expression_all(find(ismember(genes,refModel.genes)),grpIND);

disp("Getting expression based reaction scores for the reference model...needs Cobra toolbox")
expression = struct; expression.gene = genes; expression.rawValue = expression_all; expression.value = mean(expression_all,2);
[expressionRxns_GIMME, genes_rxns_asso_GIMME] = mapExpressionToReactions(refModel,expression); %% Generate reaction scores %%%

disp('Setting solver to ibm_cplex... needs IBM Cplex solver')
    setRavenSolver('ibm_cplex'); %% look RAVEN documentation to execute this !
       
disp("Method:1 :: Constraining Reference Model with GIMME at a defined threshold (twice gene expression in abs scale)");
   tic
     options=struct();
     options.solver='GIMME';
     options.threshold = rxn_cutoff_GIMME; %% Threshold set to 1; Atleast 2 times expressed in the absolute scale
     options.expressionRxns = expressionRxns_GIMME;
     options.obj_frac = 0.9;
     
     %% Partho: Change here for model reductions in terms of fluxes %%
     model_GIMME = createTissueSpecificModel(refModel,options,1); %% Change here to build a functional model; Mandatory and reduce with blocked reactions
     
     model_GIMME.id = 'HMR_2.00_GIMME';
   toc
   disp('GIMME model')
   
disp('done!')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ####### E-Flux ###### %%%%%%%%%%%%%%%%%
             
disp("Method:2 :: Next constrain the GIMME model with E-Flux algorithm. Note the model reaction's lbs and ubs are constrained with expression data... needs Cobra toolbox");
disp("Subjecting GIMME model to E-flux");

% G1 = model_GIMME.genes;
% E1 = zeros(length(model_GIMME.genes), size(grpIND,2));
% 
% for uu = 1:length(G1)
%     uu
%     idpp = find(ismember(genes,G1(uu)));
%     if(idpp>0)
%      E1(uu,:) = expression_all(idpp(1),grpIND);
%     end
%     
% end

%E1(find(ismember(model_GIMME.genes,genes)),grpIND) = expression_all(find(ismember(model_GIMME.genes,genes)),grpIND);
expression = struct; expression.gene = genes; expression.rawValue = expression_all; expression.value = mean(expression_all,2);

%expression = struct; expression.gene = G1; expression.rawValue = E1; expression.value = mean(E1,2);

[expressionRxns_EFLUX,genes_rxns_asso_Eflux] = mapExpressionToReactions(model_GIMME,expression); %% Generate reaction scores for the Healthy Obese/normal %%%

   tic
     disp('getting the reaction scores using GPR rules of the GIMME model just created....')
     expression=struct();
     expression.target = model_GIMME.rxns;
     expression.value = expressionRxns_EFLUX;
     expression.preprocessed = 1; %% constrined by reactions only!
     
     disp('Applying E-flux on GIMME constrained Model ... needs Cobra')
     CondiModel = applyEFluxConstraints(model_GIMME,expression);
     CondiModel.modelID = nom;
     CondiModel.description = strcat(nom,'_GIMME','_E-Flux');
     CondiModel.modelName = strcat(nom,'_HMR2.0_derived');
     CondiModel.id = strcat(nom,'_HMR2.0');
     
   toc
   
   disp('done!')
   disp(strcat('Created',nom));
   disp('Note this model is constrained by both first by GIMME and then by E-Fluxes....')

end