#list antibotics and target df
import pandas as pd

to_drop = ["subject_id","micro_specimen_id"]
indexnames = ['MALE', 'abnormalWCC', 'admission_location_AMBULATORY SURGERY TRANSFER', 'admission_location_CLINIC REFERRAL', 'admission_location_EMERGENCY ROOM', 'admission_location_INFORMATION NOT AVAILABLE', 'admission_location_INTERNAL TRANSFER TO OR FROM PSYCH', 'admission_location_OUTPATIENT', 'admission_location_PACU', 'admission_location_PHYSICIAN REFERRAL', 'admission_location_PROCEDURE SITE', 'admission_location_TRANSFER FROM HOSPITAL', 'admission_location_TRANSFER FROM SKILLED NURSING FACILITY', 'admission_location_WALK-IN/SELF REFERRAL', 'curr_service_CMED', 'curr_service_CSURG', 'curr_service_ENT', 'curr_service_GU', 'curr_service_GYN', 'curr_service_MED', 'curr_service_NMED', 'curr_service_NSURG', 'curr_service_OBS', 'curr_service_OMED', 'curr_service_ORTHO', 'curr_service_PSURG', 'curr_service_PSYCH', 'curr_service_SURG', 'curr_service_TRAUM', 'curr_service_TSURG', 'curr_service_UNKNOWN', 'curr_service_VSURG', 'd7AMCrx', 'd7AMKrx', 'd7AMPrx', 'd7AMXrx', 'd7ATMrx', 'd7AZMrx', 'd7CAZrx', 'd7CIPrx', 'd7CLIrx', 'd7CLRrx', 'd7CROrx', 'd7CZOrx', 'd7DAPrx', 'd7DOXrx', 'd7ERYrx', 'd7ETPrx', 'd7FEPrx', 'd7GENrx', 'd7LNZrx', 'd7MEMrx', 'd7MTRrx', 'd7NITrx', 'd7RIFrx', 'd7SAMrx', 'd7SXTrx', 'd7TOBrx', 'd7TZPrx', 'd7VANrx', 'highCRP', 'insurance_Medicaid', 'insurance_Medicare', 'insurance_Other', 'insurance_UNKNOWN', 'marital_status_DIVORCED', 'marital_status_MARRIED', 'marital_status_SINGLE', 'marital_status_UNKNOWN', 'marital_status_WIDOWED', 'ob_freq', 'pAMCrx', 'pAMKi', 'pAMKnt', 'pAMKr', 'pAMKrx', 'pAMKs', 'pAMPc', 'pAMPcS', 'pAMPi', 'pAMPnt', 'pAMPr', 'pAMPrx', 'pAMPs', 'pAMXrx', 'pATMrx', 'pAZMrx', 'pCATH', 'pCAZi', 'pCAZnt', 'pCAZr', 'pCAZrx', 'pCAZs', 'pCIPi', 'pCIPnt', 'pCIPr', 'pCIPrx', 'pCIPs', 'pCLIi', 'pCLInt', 'pCLIr', 'pCLIrx', 'pCLIs', 'pCLRrx', 'pCROi', 'pCROnt', 'pCROr', 'pCROrx', 'pCROs', 'pCZOi', 'pCZOnt', 'pCZOr', 'pCZOrx', 'pCZOs', 'pChemo', 'pDAPrx', 'pDISC', 'pDNR', 'pDOXrx', 'pERYr', 'pERYrx', 'pETPrx', 'pFEPi', 'pFEPnt', 'pFEPr', 'pFEPrx', 'pFEPs', 'pGENi', 'pGENnt', 'pGENr', 'pGENrx', 'pGENs', 'pHADM', 'pHyd', 'pICD_A', 'pICD_B', 'pICD_C', 'pICD_D', 'pICD_E', 'pICD_F', 'pICD_G', 'pICD_H', 'pICD_I', 'pICD_J', 'pICD_K', 'pICD_L', 'pICD_M', 'pICD_N', 'pICD_O', 'pICD_P', 'pICD_Q', 'pICD_R', 'pICD_S', 'pICD_T', 'pICD_U', 'pICD_V', 'pICD_W', 'pICD_X', 'pICD_Y', 'pICD_Z', 'pICU', 'pLNZrx', 'pLVXi', 'pLVXnt', 'pLVXr', 'pLVXs', 'pMEMi', 'pMEMnt', 'pMEMr', 'pMEMrx', 'pMEMs', 'pMTRrx', 'pNGT', 'pNH', 'pNITi', 'pNITnt', 'pNITr', 'pNITrx', 'pNITs', 'pNUTR', 'pNeph', 'pOT', 'pObese', 'pOverweight', 'pPENi', 'pPENnt', 'pPENr', 'pPENs', 'pPROC_0', 'pPROC_1', 'pPROC_2', 'pPROC_3', 'pPROC_4', 'pPROC_5', 'pPROC_6', 'pPROC_7', 'pPROC_8', 'pPROC_9', 'pPROC_A', 'pPROC_B', 'pPROC_C', 'pPROC_D', 'pPROC_E', 'pPROC_F', 'pPROC_G', 'pPROC_H', 'pPROC_I', 'pPROC_J', 'pPROC_K', 'pPROC_L', 'pPROC_M', 'pPROC_N', 'pPROC_O', 'pPROC_P', 'pPROC_Q', 'pPROC_R', 'pPROC_S', 'pPROC_T', 'pPROC_X', 'pPhysio', 'pPsych', 'pRIFrx', 'pRestr', 'pSAMi', 'pSAMnt', 'pSAMr', 'pSAMrx', 'pSAMs', 'pSXTi', 'pSXTnt', 'pSXTr', 'pSXTrx', 'pSXTs', 'pSocial', 'pSurg', 'pTCYi', 'pTCYnt', 'pTCYr', 'pTCYs', 'pTOBi', 'pTOBnt', 'pTOBr', 'pTOBrx', 'pTOBs', 'pTPN', 'pTZPi', 'pTZPnt', 'pTZPr', 'pTZPrx', 'pTZPs', 'pUnderweight', 'pVANi', 'pVANnt', 'pVANr', 'pVANrx', 'pVANs', 'provider_id', 'race_AMERICAN INDIAN/ALASKA NATIVE', 'race_ASIAN', 'race_ASIAN - ASIAN INDIAN', 'race_ASIAN - CHINESE', 'race_ASIAN - KOREAN', 'race_ASIAN - SOUTH EAST ASIAN', 'race_BLACK/AFRICAN', 'race_BLACK/AFRICAN AMERICAN', 'race_BLACK/CAPE VERDEAN', 'race_BLACK/CARIBBEAN ISLAND', 'race_HISPANIC OR LATINO', 'race_HISPANIC/LATINO - CENTRAL AMERICAN', 'race_HISPANIC/LATINO - COLUMBIAN', 'race_HISPANIC/LATINO - CUBAN', 'race_HISPANIC/LATINO - DOMINICAN', 'race_HISPANIC/LATINO - GUATEMALAN', 'race_HISPANIC/LATINO - HONDURAN', 'race_HISPANIC/LATINO - MEXICAN', 'race_HISPANIC/LATINO - PUERTO RICAN', 'race_HISPANIC/LATINO - SALVADORAN', 'race_MULTIPLE RACE/ETHNICITY', 'race_NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER', 'race_OTHER', 'race_PATIENT DECLINED TO ANSWER', 'race_PORTUGUESE', 'race_SOUTH AMERICAN', 'race_UNABLE TO OBTAIN', 'race_UNKNOWN', 'race_WHITE', 'race_WHITE - BRAZILIAN', 'race_WHITE - EASTERN EUROPEAN', 'race_WHITE - OTHER EUROPEAN', 'race_WHITE - RUSSIAN', 'standard_age']
target_df = "daily_urines.csv"
probs_cols = ['Antimicrobial', 'I', 'R', 'S', 'Variable_tuning', 'micro_specimen_id', 'subject_id']
ref_df = pd.read_csv("urines5.csv")
ref_df['standard_age'] = ref_df['standard_age'].map(str)
indexnames = sorted(pd.get_dummies(ref_df))
#Create empty results dataframe
probs_df_overall = pd.DataFrame()
probs_df_overall = probs_df_overall.reindex(columns=['Antimicrobial', 'I', 'R', 'S','NT', 'Variable_tuning', 'micro_specimen_id', 'subject_id'])

#####################################


#AMPICILLIN

#####################################
prob_predict('Ampicillin','AMP',to_drop,target_df,'ampfile.pickle','ampfits.pickle',
             full_varlist=indexnames,
             result_optimised_for="R",ref=ref_df)


probs_df_final = probs_df_final.reindex(columns=['Antimicrobial', 'I', 'R', 'S','NT', 'Variable_tuning', 'micro_specimen_id', 'subject_id'])
probs_df_overall = pd.concat([probs_df_overall,probs_df_final])


#####################################

#AMPICILLIN-SULBACTAM

#####################################

prob_predict('Ampicillin-sulbactam','SAM',to_drop,target_df,'samfile.pickle','samfits.pickle',
             full_varlist=indexnames,
             result_optimised_for="R",ref=ref_df)

probs_df_final = probs_df_final.reindex(columns=['Antimicrobial', 'I', 'R', 'S','NT', 'Variable_tuning', 'micro_specimen_id', 'subject_id'])
probs_df_overall = pd.concat([probs_df_overall,probs_df_final])


#####################################

#PIPERACILLIN-TAZOBACTAM

#####################################

prob_predict('Piperacillin-tazobactam','TZP',to_drop,target_df,'tzpfile.pickle','tzpfits.pickle',
             full_varlist=indexnames,
             result_optimised_for="R",ref=ref_df)

probs_df_final = probs_df_final.reindex(columns=['Antimicrobial', 'I', 'R', 'S','NT', 'Variable_tuning', 'micro_specimen_id', 'subject_id'])
probs_df_overall = pd.concat([probs_df_overall,probs_df_final])


#####################################

#CEFAZOLIN

#####################################

prob_predict('Cefazolin','CZO',to_drop,target_df,'czofile.pickle','czofits.pickle',
             full_varlist=indexnames,
             result_optimised_for="R",ref=ref_df)

probs_df_final = probs_df_final.reindex(columns=['Antimicrobial', 'I', 'R', 'S','NT', 'Variable_tuning', 'micro_specimen_id', 'subject_id'])
probs_df_overall = pd.concat([probs_df_overall,probs_df_final])

#####################################

#CEFTRIAXONE

#####################################

prob_predict('Ceftriaxone','CRO',to_drop,target_df,'crofile.pickle','crofits.pickle',
             full_varlist=indexnames,
             result_optimised_for="R",ref=ref_df)

probs_df_final = probs_df_final.reindex(columns=['Antimicrobial', 'I', 'R', 'S','NT', 'Variable_tuning', 'micro_specimen_id', 'subject_id'])
probs_df_overall = pd.concat([probs_df_overall,probs_df_final])

#####################################

#CEFTAZIDIME

#####################################

prob_predict('Ceftazidime','CAZ',to_drop,target_df,'cazfile.pickle','cazfits.pickle',
             full_varlist=indexnames,
             result_optimised_for="R",ref=ref_df)

probs_df_final = probs_df_final.reindex(columns=['Antimicrobial', 'I', 'R', 'S','NT', 'Variable_tuning', 'micro_specimen_id', 'subject_id'])
probs_df_overall = pd.concat([probs_df_overall,probs_df_final])

#####################################

#CEFEPIME

#####################################

prob_predict('Cefepime','FEP',to_drop,target_df,'fepfile.pickle','fepfits.pickle',
             full_varlist=indexnames,
             result_optimised_for="R",ref=ref_df)

probs_df_final = probs_df_final.reindex(columns=['Antimicrobial', 'I', 'R', 'S','NT', 'Variable_tuning', 'micro_specimen_id', 'subject_id'])
probs_df_overall = pd.concat([probs_df_overall,probs_df_final])

#####################################

#MEROPENEM

#####################################

prob_predict('Meropenem','MEM',to_drop,target_df,'memfile.pickle','memfits.pickle',
             full_varlist=indexnames,
             result_optimised_for="R",ref=ref_df)

probs_df_final = probs_df_final.reindex(columns=['Antimicrobial', 'I', 'R', 'S','NT', 'Variable_tuning', 'micro_specimen_id', 'subject_id'])
probs_df_overall = pd.concat([probs_df_overall,probs_df_final])

#####################################

#CIPROFLOXACIN

#####################################

prob_predict('Ciprofloxacin','CIP',to_drop,target_df,'cipfile.pickle','cipfits.pickle',
             full_varlist=indexnames,
             result_optimised_for="R",ref=ref_df)

probs_df_final = probs_df_final.reindex(columns=['Antimicrobial', 'I', 'R', 'S','NT', 'Variable_tuning', 'micro_specimen_id', 'subject_id'])
probs_df_overall = pd.concat([probs_df_overall,probs_df_final])

#####################################

#GENTAMICIN

#####################################

prob_predict('Gentamicin','GEN',to_drop,target_df,'genfile.pickle','genfits.pickle',
             full_varlist=indexnames,
             result_optimised_for="R",ref=ref_df)

probs_df_final = probs_df_final.reindex(columns=['Antimicrobial', 'I', 'R', 'S','NT', 'Variable_tuning', 'micro_specimen_id', 'subject_id'])
probs_df_overall = pd.concat([probs_df_overall,probs_df_final])

#####################################

#TRIMETHOPRIM-SULFAMETHOXAZOLE

#####################################

prob_predict('Trimethoprim-sulfamethoxazole','SXT',to_drop,target_df,'sxtfile.pickle','sxtfits.pickle',
             full_varlist=indexnames,
             result_optimised_for="R",ref=ref_df)

probs_df_final = probs_df_final.reindex(columns=['Antimicrobial', 'I', 'R', 'S','NT', 'Variable_tuning', 'micro_specimen_id', 'subject_id'])
probs_df_overall = pd.concat([probs_df_overall,probs_df_final])

#####################################

#NITROFURANTOIN

#####################################

prob_predict('Nitrofurantoin','NIT',to_drop,target_df,'nitfile.pickle','nitfits.pickle',
             full_varlist=indexnames,
             result_optimised_for="R",ref=ref_df)

probs_df_final = probs_df_final.reindex(columns=['Antimicrobial', 'I', 'R', 'S','NT', 'Variable_tuning', 'micro_specimen_id', 'subject_id'])
probs_df_overall = pd.concat([probs_df_overall,probs_df_final])

#####################################


#WRITE TO CSV
probs_df_overall.to_csv("probs_df_overall.csv")





