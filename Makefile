## ### Make the rescuability files ## ###
InputDIR = ./File/Original_SData/
OutputDIR = ./File/RescuabilityFile/

python:
            python ./MakeRescuabilityFile.py InputDIR OutputDIR
            
## ### Make the conservstion scores for each site files ## ###
InputDIR = ./File/Original_SData/
OutputDIR = ./File/RSAWCNConservationInterfaceFile/

python:
            python ./MakeSingleSiteConservationFile.py InputDIR OutputDIR
            
## ### Make the rescuability score for each site files ## ###
InputDIR_1 = ./File/RescuabilityFile/
InputDIR_2 = ./File/RSAWCNConservationInterfaceFile/
OutputDIR = $(InputDIR_1)

python:
        python ./MakeSingleSiteRescuabilityFile.py InputDIR_1 InputDIR_2 OutputDIR

## ### Make the ML files ## ###
InputDIR_1 = ./File/RescuabilityFile/
InputDIR_2 = ./File/RSAWCNConservationInterfaceFile/
InputDIR_3 = ./File/AAindex/
OutputDIR = ./File/MLFile/

python:
        python ./MakeMachineLearningFile.py InputDIR_1 InputDIR_2 InputDIR_3 OutputDIR

## ### Make the features predict fitness files ## ###
InputDIR_1 = ./File/AAindex/
InputDIR_2 = ./File/Original_SData/
OutputDIR = ./File/FeaturePredictFitness/

python:
        python ./MakeFeaturePredictFitnessFile.py InputDIR_1 InputDIR_2 OutputDIR

## ### Make the features predict compensation ## ###
InputDIR_1 = ./File/FeaturePredictFitness/
InputDIR_2 = ./File/RescuabilityFile/
InputDIR_3 = ./File/Original_SData/
InputDIR_4 = ./File/AAindex/
OutputDIR = ./File/FeaturePredictCompensation/

python:
        python ./MakeFeaturePredictCompensationSuperset.py InputDIR_1 InputDIR_2 InputDIR_3 InputDIR_4 OutputDIR
   
## ### Make the features predict compensation for subset files ## ###
InputDIR_1 = ./File/FeaturePredictFitness/
InputDIR_2 = ./File/RescuabilityFile/
InputDIR_3 = ./File/Original_SData/
InputDIR_4 = ./File/AAindex/
InputDIR_5 = ./File/FeaturePredictCompensation/
OutputDIR = ${InputDIR_5}

python:
        python ./MakeFeaturePredictCompensationSubset.py InputDIR_1 InputDIR_2 InputDIR_3 InputDIR_4 InputDIR_5 OutputDIR

## ### Make the SC files ## ###
InputDIR_1 = ./File/RescuabilityFile/
InputDIR_2 = ./File/Original_SData/
OutputDIR = ./File/SuperCompensationFile/

python:
        python ./MakeSubsCompensateAlleles.py InputDIR_1 InputDIR_2  OutputDIR

## ### Make the buffering ability of SC files ## ###
InputDIR = ./File/SuperCompensationFile/
OutputDIR = ${InputDIR}

python:
        python ./MakeSCBufferingAbility.py InputDIR OutputDIR
        
echo "Thanks for using these scripts!" 

        





