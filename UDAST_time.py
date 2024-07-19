#load in, check and clean

import pandas as pd

drops = ["AMP","SAM","TZP","CZO","CRO","CAZ","FEP","MEM","CIP","GEN","SXT","NIT","VAN" ]
ref_df = pd.read_csv("urines_5t_1.csv")
ref_df['standard_age'] = ref_df['standard_age'].map(str)
dummy_df = pd.get_dummies(ref_df.drop(drops, axis=1))
vari_overall = {}
c_values = {}
to_drop = ["subject_id","micro_specimen_id"]
class_reps = {}
prot_vars = []
gen_cols = ['MALE']
age_cols = [col for col in dummy_df.columns if 'standard_age' in col]
race_cols = [col for col in dummy_df.columns if 'race' in col]
language_cols = [col for col in dummy_df.columns if 'language' in col]
insurance_cols = [col for col in dummy_df.columns if 'insurance' in col]
marital_cols = [col for col in dummy_df.columns if 'marital' in col]
prot_vars.extend(gen_cols)
prot_vars.extend(age_cols)
prot_vars.extend(race_cols)
prot_vars.extend(language_cols)
prot_vars.extend(insurance_cols)
prot_vars.extend(marital_cols)


#####WITHIN OWN TIME PERIOD

def own_time_per(filename_t):

    global aucrocs
    aucrocs = {}

    for i in list(range(1, 21, 1)):


        rocs = {}

        ###############################

        #AMPICILLIN

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['AMP']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'AMP', y)


        #set Y variable
        urines5['Y'] = urines5['AMP']
        urines5 = urines5.drop('AMP',axis=1)
        features = urines5.drop('Y',axis=1)


        #C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                    target=urines5["Y"],
                    ht_features=features,
                    test=0.2,
                    random=i,
                    folds=6,
                    model_type="ovr",
                    cw='balanced',
                    scorer='roc_auc_ovr',
                    target_ab = 'Ampicillin')

        vari_overall['AMP'] = vari_list
        c_values['AMP'] = c_value

        #LR final run
        LR_multi_final(target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=i,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Ampicillin",
                 av="micro")

        class_reps['AMP'] = class_report

        #save fit, variables, and hyperparameters
        filename = "amp_" + filename_t + ".pickle"
        fitname = "amp_" +filename_t + ".pickle"
        with open(filename,'wb') as f:
            pickle.dump(vari_list,f)
        with open(fitname,'wb') as f:
            pickle.dump(model_dict,f)

        rocs['AMP'] = roc_auc


        ###############################

        #AMPICILLIN-SULBACTAM

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['SAM']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'SAM', y)

        #set Y variable
        urines5['Y'] = urines5['SAM']
        urines5 = urines5.drop('SAM',axis=1)
        features = urines5.drop('Y',axis=1)

        #C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                    target=urines5["Y"],
                    ht_features=features,
                    test=0.2,
                    random=i,
                    folds=6,
                    model_type="ovr",
                    cw='balanced',
                    scorer='roc_auc_ovr',
                    target_ab = 'Ampicillin-sulbactam')

        vari_overall['SAM'] = vari_list
        c_values['SAM'] = c_value

        #LR final run
        LR_multi_final(target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=i,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Ampicillin-sulbactam",
                 av="micro")

        class_reps['SAM'] = class_report

        #save fit, variables, and hyperparameters
        filename = "sam_" + filename_t + ".pickle"
        fitname = "sam_" + filename_t + ".pickle"
        with open(filename,'wb') as f:
            pickle.dump(vari_list,f)
        with open(fitname,'wb') as f:
            pickle.dump(model_dict,f)

        rocs['SAM'] = roc_auc

        ###############################

        #PIPERACILLIN-TAZOBACTAM

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['TZP']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'TZP', y)

        #set Y variable
        urines5['Y'] = urines5['TZP']
        urines5 = urines5.drop('TZP',axis=1)
        features = urines5.drop('Y',axis=1)

        #C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                    target=urines5["Y"],
                    ht_features=features,
                    test=0.2,
                    random=i,
                    folds=6,
                    model_type="ovr",
                    cw='balanced',
                    scorer='roc_auc_ovr',
                    target_ab = 'Piperacillin-tazobactam')

        vari_overall['TZP'] = vari_list
        c_values['TZP'] = c_value

        #LR final run
        LR_multi_final(target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=i,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Piperacillin-tazobactam",
                 av="micro")

        class_reps['TZP'] = class_report

        #save fit, variables, and hyperparameters
        filename = "tzp_" + filename_t + ".pickle"
        fitname = "tzp_" + filename_t + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['TZP'] = roc_auc

        ###############################

        #CEFAZOLIN

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['CZO']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'CZO', y)

        #set Y variable
        urines5['Y'] = urines5['CZO']
        urines5 = urines5.drop('CZO',axis=1)
        features = urines5.drop('Y',axis=1)

        #C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                    target=urines5["Y"],
                    ht_features=features,
                    test=0.2,
                    random=i,
                    folds=6,
                    model_type="ovr",
                    cw='balanced',
                    scorer='roc_auc_ovr',
                    target_ab = 'Cefazolin')

        vari_overall['CZO'] = vari_list
        c_values['CZO'] = c_value

        #LR final run
        LR_multi_final(target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=i,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Cefazolin",
                 av="micro")

        class_reps['CZO'] = class_report

        #save fit, variables, and hyperparameters
        filename = "czo_" + filename_t + ".pickle"
        fitname = "czo_" + filename_t + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['CZO'] = roc_auc

        ###############################

        #CEFTRIAXONE

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['CRO']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'CRO', y)

        #set Y variable
        urines5['Y'] = urines5['CRO']
        urines5 = urines5.drop('CRO',axis=1)
        features = urines5.drop('Y',axis=1)

        #C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                    target=urines5["Y"],
                    ht_features=features,
                    test=0.2,
                    random=i,
                    folds=6,
                    model_type="ovr",
                    cw='balanced',
                    scorer='roc_auc_ovr',
                    target_ab='Ceftriaxone')

        vari_overall['CRO'] = vari_list
        c_values['CRO'] = c_value

        #LR final run
        LR_multi_final(target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=i,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Ceftriaxone",
                 av="micro")

        class_reps['CRO'] = class_report

        #save fit, variables, and hyperparameters
        filename = "cro_" + filename_t + ".pickle"
        fitname = "cro_" + filename_t + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['CRO'] = roc_auc

        ###############################

        #CEFTAZIDIME

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['CAZ']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'CAZ', y)

        #set Y variable
        urines5['Y'] = urines5['CAZ']
        urines5 = urines5.drop('CAZ',axis=1)
        features = urines5.drop('Y',axis=1)

        #C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                    target=urines5["Y"],
                    ht_features=features,
                    test=0.2,
                    random=i,
                    folds=6,
                    model_type="ovr",
                    cw='balanced',
                    scorer='roc_auc_ovr',
                    target_ab = 'Ceftazidime')

        vari_overall['CAZ'] = vari_list
        c_values['CAZ'] = c_value

        #LR final run
        LR_multi_final(target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=i,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Ceftazidime",
                 av="micro")


        class_reps['CAZ'] = class_report

        #save fit, variables, and hyperparameters
        filename = "caz_" + filename_t + ".pickle"
        fitname = "caz_" + filename_t + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['CAZ'] = roc_auc

        ###############################

        #CEFEPIME

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['FEP']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'FEP', y)

        #set Y variable
        urines5['Y'] = urines5['FEP']
        urines5 = urines5.drop('FEP',axis=1)
        features = urines5.drop('Y',axis=1)

        #C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                    target=urines5["Y"],
                    ht_features=features,
                    test=0.2,
                    random=i,
                    folds=6,
                    model_type="ovr",
                    cw='balanced',
                    scorer='roc_auc_ovr',
                    target_ab = 'Cefepime')

        vari_overall['FEP'] = vari_list
        c_values['FEP'] = c_value

        #LR final run
        LR_multi_final(target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=i,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Cefepime",
                 av="micro")

        class_reps['FEP'] = class_report

        #save fit, variables, and hyperparameters
        filename = "fep_" + filename_t + ".pickle"
        fitname = "fep_" + filename_t + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['FEP'] = roc_auc

        ###############################

        #MEROPENEM

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['MEM']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'MEM', y)

        #set Y variable
        urines5['Y'] = urines5['MEM']
        urines5 = urines5.drop('MEM',axis=1)
        features = urines5.drop('Y',axis=1)

        #C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                    target=urines5["Y"],
                    ht_features=features,
                    test=0.2,
                    random=i,
                    folds=6,
                    model_type="ovr",
                    cw='balanced',
                    scorer='roc_auc_ovr',
                    target_ab = 'Meropenem')

        vari_overall['MEM'] = vari_list
        c_values['MEM'] = c_value

        #LR final run
        LR_multi_final(target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=i,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Meropenem",
                 av="micro")

        class_reps['MEM'] = class_report

        #save fit, variables, and hyperparameters
        filename = "mem_" + filename_t + ".pickle"
        fitname = "mem_" + filename_t + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['MEM'] = roc_auc

        ###############################

        #CIPROFLOXACIN

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['CIP']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'CIP', y)

        #set Y variable
        urines5['Y'] = urines5['CIP']
        urines5 = urines5.drop('CIP',axis=1)
        features = urines5.drop('Y',axis=1)

        #C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                    target=urines5["Y"],
                    ht_features=features,
                    test=0.2,
                    random=i,
                    folds=6,
                    model_type="ovr",
                    cw='balanced',
                    scorer='roc_auc_ovr',
                    target_ab = 'Ciprofloxacin')

        vari_overall['CIP'] = vari_list
        c_values['CIP'] = c_value

        #LR final run
        LR_multi_final(target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=i,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Ciprofloxacin",
                 av="micro")

        class_reps['CIP'] = class_report

        #save fit, variables, and hyperparameters
        filename = "cip_" + filename_t + ".pickle"
        fitname = "cip_" + filename_t + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['CIP'] = roc_auc

        ###############################

        #GENTAMICIN

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['GEN']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'GEN', y)

        #set Y variable
        urines5['Y'] = urines5['GEN']
        urines5 = urines5.drop('GEN',axis=1)
        features = urines5.drop('Y',axis=1)

        #C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                    target=urines5["Y"],
                    ht_features=features,
                    test=0.2,
                    random=i,
                    folds=6,
                    model_type="ovr",
                    cw='balanced',
                    scorer='roc_auc_ovr',
                    target_ab = 'Gentamicin')

        vari_overall['GEN'] = vari_list
        c_values['GEN'] = c_value

        #LR final run
        LR_multi_final(target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=i,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Gentamicin",
                 av="micro")

        class_reps['GEN'] = class_report

        #save fit, variables, and hyperparameters
        filename = "gen_" + filename_t + ".pickle"
        fitname = "gen_" + filename_t + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['GEN'] = roc_auc

        ###############################

        #TRIMETHOPRIM-SULFAMETHOXAZOLE

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['SXT']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'SXT', y)

        #set Y variable
        urines5['Y'] = urines5['SXT']
        urines5 = urines5.drop('SXT',axis=1)
        features = urines5.drop('Y',axis=1)

        #C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                    target=urines5["Y"],
                    ht_features=features,
                    test=0.2,
                    random=i,
                    folds=6,
                    model_type="ovr",
                    cw='balanced',
                    scorer='roc_auc_ovr',
                    target_ab = 'Trimethoprim-sulfamethoxazole',
                    targ_result = 'R')

        vari_overall['SXT'] = vari_list
        c_values['SXT'] = c_value

        #LR final run
        LR_multi_final(target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=i,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Trimethoprim-sulfamethoxazole",
                 av="micro")

        class_reps['SXT'] = class_report

        #save fit, variables, and hyperparameters
        filename = "sxt_" + filename_t + ".pickle"
        fitname = "sxt_" + filename_t + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['SXT'] = roc_auc

        ###############################

        #NITROFURANTOIN

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['NIT']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'NIT', y)

        #set Y variable
        urines5['Y'] = urines5['NIT']
        urines5 = urines5.drop('NIT',axis=1)
        features = urines5.drop('Y',axis=1)

        #C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                    target=urines5["Y"],
                    ht_features=features,
                    test=0.2,
                    random=i,
                    folds=6,
                    model_type="ovr",
                    cw='balanced',
                    scorer='roc_auc_ovr',
                    target_ab = 'Nitrofurantoin')

        vari_overall['NIT'] = vari_list
        c_values['NIT'] = c_value

        #LR final run
        LR_multi_final(target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=i,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Nitrofurantoin",
                 av="micro")

        class_reps['NIT'] = class_report

        #save fit, variables, and hyperparameters
        filename = "nit_" + filename_t + ".pickle"
        fitname = "nit_" + filename_t + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['NIT'] = roc_auc

        ###############################

        #VANCOMYCIN

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['VAN']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'VAN', y)

        #set Y variable
        urines5['Y'] = urines5['VAN']
        urines5 = urines5.drop('VAN',axis=1)
        features = urines5.drop('Y',axis=1)

        #C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                    target=urines5["Y"],
                    ht_features=features,
                    test=0.2,
                    random=i,
                    folds=6,
                    model_type="ovr",
                    cw='balanced',
                    scorer='roc_auc_ovr',
                    target_ab = 'Vancomycin',
                    targ_result='R')

        vari_overall['VAN'] = vari_list
        c_values['VAN'] = c_value

        #LR final run
        LR_multi_final(target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=i,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Vancomycin",
                 av="micro")

        class_reps['VAN'] = class_report

        #save fit, variables, and hyperparameters
        filename = "van_" + filename_t + ".pickle"
        fitname = "van_" + filename_t + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['VAN'] = roc_auc

        aucrocs[i] = rocs

own_time_per("urines_5t_1.csv")
aucrocs_5t_1_1 = aucrocs
aucrocs_5t_1_1['df1'] = 1
aucrocs_5t_1_1['df2'] = 1
aucrocs_5t_1_1 = pd.DataFrame(aucrocs_5t_1_1)
own_time_per("urines_5t_2.csv")
aucrocs_5t_2_2 = aucrocs
aucrocs_5t_2_2['df1'] = 2
aucrocs_5t_2_2['df2'] = 2
aucrocs_5t_2_2 = pd.DataFrame(aucrocs_5t_2_2)
own_time_per("urines_5t_3.csv")
aucrocs_5t_3_3 = aucrocs
aucrocs_5t_3_3['df1'] = 3
aucrocs_5t_3_3['df2'] = 3
aucrocs_5t_3_3 = pd.DataFrame(aucrocs_5t_3_3)
own_time_per("urines_5t_4.csv")
aucrocs_5t_4_4 = aucrocs
aucrocs_5t_4_4['df1'] = 4
aucrocs_5t_4_4['df2'] = 4
aucrocs_5t_4_4 = pd.DataFrame(aucrocs_5t_4_4)




#######ACROSS TIME PERIODS

def across_time_per(filename_t,filename_t2):

    global aucrocs
    aucrocs={}

    for i in list(range(1,21,1)):

        rocs = {}

        ###############################

        # AMPICILLIN

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['AMP']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'AMP', y)

        # set Y variable
        urines5['Y'] = urines5['AMP']
        urines5 = urines5.drop('AMP', axis=1)
        features = urines5.drop('Y', axis=1)

        urines6 = pd.read_csv(filename_t2)
        urines6['standard_age'] = urines6['standard_age'].map(str)
        y = urines6['AMP']
        urines6 = pd.get_dummies(urines6.drop(drops, axis=1))
        urines6.insert(0, 'AMP', y)

        # set Y variable
        urines6['Y'] = urines6['AMP']
        urines6 = urines6.drop('AMP', axis=1)
        features2 = urines6.drop('Y', axis=1)

        missing_columns = set(features.columns) - set(features2.columns)
        for col in missing_columns:
            features2[col] = False

        features2 = features2[features.columns]

        # C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                          target=urines5["Y"],
                          ht_features=features,
                          test=0.2,
                          random=i,
                          folds=6,
                          model_type="ovr",
                          cw='balanced',
                          scorer='roc_auc_ovr',
                          target_ab='Ampicillin')

        vari_overall['AMP'] = vari_list
        c_values['AMP'] = c_value

        # LR final run
        LR_multi_final_t(target=urines5["Y"],
                       final_features=features,
                       test=0.2,
                       random=i,
                       C_val=c_value,
                       cw='balanced',
                       model_type="ovr",
                       thr=0.5,
                       ab_of_interest="Ampicillin",
                       av="micro",
                       target2=urines6["Y"],
                       test_datf=features2)

        class_reps['AMP'] = class_report

        # save fit, variables, and hyperparameters
        filename = "amp_" + filename_t + "_" + filename_t2 + ".pickle"
        fitname = "amp_" + filename_t + "_" + filename_t2 + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['AMP'] = roc_auc

        ###############################

        # AMPICILLIN-SULBACTAM

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['SAM']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'SAM', y)

        # set Y variable
        urines5['Y'] = urines5['SAM']
        urines5 = urines5.drop('SAM', axis=1)
        features = urines5.drop('Y', axis=1)

        urines6 = pd.read_csv(filename_t2)
        urines6['standard_age'] = urines6['standard_age'].map(str)
        y = urines6['SAM']
        urines6 = pd.get_dummies(urines6.drop(drops, axis=1))
        urines6.insert(0, 'SAM', y)

        # set Y variable
        urines6['Y'] = urines6['SAM']
        urines6 = urines6.drop('SAM', axis=1)
        features2 = urines6.drop('Y', axis=1)

        missing_columns = set(features.columns) - set(features2.columns)
        for col in missing_columns:
            features2[col] = False

        features2 = features2[features.columns]

        # C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                          target=urines5["Y"],
                          ht_features=features,
                          test=0.2,
                          random=i,
                          folds=6,
                          model_type="ovr",
                          cw='balanced',
                          scorer='roc_auc_ovr',
                          target_ab='Ampicillin-sulbactam')

        vari_overall['SAM'] = vari_list
        c_values['SAM'] = c_value

        # LR final run
        LR_multi_final_t(target=urines5["Y"],
                       final_features=features,
                       test=0.2,
                       random=i,
                       C_val=c_value,
                       cw='balanced',
                       model_type="ovr",
                       thr=0.5,
                       ab_of_interest="Ampicillin-sulbactam",
                       av="micro",
                       target2=urines6["Y"],
                       test_datf=features2)

        class_reps['SAM'] = class_report

        # save fit, variables, and hyperparameters
        filename = "sam_" + filename_t + "_" + filename_t2 + ".pickle"
        fitname = "sam_" + filename_t + "_" + filename_t2 + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['SAM'] = roc_auc

        ###############################

        # PIPERACILLIN-TAZOBACTAM

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['TZP']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'TZP', y)

        # set Y variable
        urines5['Y'] = urines5['TZP']
        urines5 = urines5.drop('TZP', axis=1)
        features = urines5.drop('Y', axis=1)

        urines6 = pd.read_csv(filename_t2)
        urines6['standard_age'] = urines6['standard_age'].map(str)
        y = urines6['TZP']
        urines6 = pd.get_dummies(urines6.drop(drops, axis=1))
        urines6.insert(0, 'TZP', y)

        # set Y variable
        urines6['Y'] = urines6['TZP']
        urines6 = urines6.drop('TZP', axis=1)
        features2 = urines6.drop('Y', axis=1)

        missing_columns = set(features.columns) - set(features2.columns)
        for col in missing_columns:
            features2[col] = False

        features2 = features2[features.columns]

        # C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                          target=urines5["Y"],
                          ht_features=features,
                          test=0.2,
                          random=i,
                          folds=6,
                          model_type="ovr",
                          cw='balanced',
                          scorer='roc_auc_ovr',
                          target_ab='Piperacillin-tazobactam')

        vari_overall['TZP'] = vari_list
        c_values['TZP'] = c_value

        # LR final run
        LR_multi_final_t(target=urines5["Y"],
                       final_features=features,
                       test=0.2,
                       random=i,
                       C_val=c_value,
                       cw='balanced',
                       model_type="ovr",
                       thr=0.5,
                       ab_of_interest="Piperacillin-tazobactam",
                       av="micro",
                       target2=urines6["Y"],
                       test_datf=features2)

        class_reps['TZP'] = class_report

        # save fit, variables, and hyperparameters
        filename = "tzp_" + filename_t + "_" + filename_t2 + ".pickle"
        fitname = "tzp_" + filename_t + "_" + filename_t2 + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['TZP'] = roc_auc

        ###############################

        # CEFAZOLIN

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['CZO']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'CZO', y)

        # set Y variable
        urines5['Y'] = urines5['CZO']
        urines5 = urines5.drop('CZO', axis=1)
        features = urines5.drop('Y', axis=1)

        urines6 = pd.read_csv(filename_t2)
        urines6['standard_age'] = urines6['standard_age'].map(str)
        y = urines6['CZO']
        urines6 = pd.get_dummies(urines6.drop(drops, axis=1))
        urines6.insert(0, 'CZO', y)

        # set Y variable
        urines6['Y'] = urines6['CZO']
        urines6 = urines6.drop('CZO', axis=1)
        features2 = urines6.drop('Y', axis=1)

        missing_columns = set(features.columns) - set(features2.columns)
        for col in missing_columns:
            features2[col] = False

        features2 = features2[features.columns]

        # C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                          target=urines5["Y"],
                          ht_features=features,
                          test=0.2,
                          random=i,
                          folds=6,
                          model_type="ovr",
                          cw='balanced',
                          scorer='roc_auc_ovr',
                          target_ab='Cefazolin')

        vari_overall['CZO'] = vari_list
        c_values['CZO'] = c_value

        # LR final run
        LR_multi_final_t(target=urines5["Y"],
                       final_features=features,
                       test=0.2,
                       random=i,
                       C_val=c_value,
                       cw='balanced',
                       model_type="ovr",
                       thr=0.5,
                       ab_of_interest="Cefazolin",
                       av="micro",
                       target2=urines6["Y"],
                       test_datf=features2)

        class_reps['CZO'] = class_report

        # save fit, variables, and hyperparameters
        filename = "czo_" + filename_t + "_" + filename_t2 + ".pickle"
        fitname = "czo_" + filename_t + "_" + filename_t2 + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['CZO'] = roc_auc

        ###############################

        # CEFTRIAXONE

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['CRO']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'CRO', y)

        # set Y variable
        urines5['Y'] = urines5['CRO']
        urines5 = urines5.drop('CRO', axis=1)
        features = urines5.drop('Y', axis=1)

        urines6 = pd.read_csv(filename_t2)
        urines6['standard_age'] = urines6['standard_age'].map(str)
        y = urines6['CRO']
        urines6 = pd.get_dummies(urines6.drop(drops, axis=1))
        urines6.insert(0, 'CRO', y)

        # set Y variable
        urines6['Y'] = urines6['CRO']
        urines6 = urines6.drop('CRO', axis=1)
        features2 = urines6.drop('Y', axis=1)

        missing_columns = set(features.columns) - set(features2.columns)
        for col in missing_columns:
            features2[col] = False

        features2 = features2[features.columns]

        # C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                          target=urines5["Y"],
                          ht_features=features,
                          test=0.2,
                          random=i,
                          folds=6,
                          model_type="ovr",
                          cw='balanced',
                          scorer='roc_auc_ovr',
                          target_ab='Ceftriaxone')

        vari_overall['CRO'] = vari_list
        c_values['CRO'] = c_value

        # LR final run
        LR_multi_final_t(target=urines5["Y"],
                       final_features=features,
                       test=0.2,
                       random=i,
                       C_val=c_value,
                       cw='balanced',
                       model_type="ovr",
                       thr=0.5,
                       ab_of_interest="Ceftriaxone",
                       av="micro",
                       target2=urines6["Y"],
                       test_datf=features2)

        class_reps['CRO'] = class_report

        # save fit, variables, and hyperparameters
        filename = "cro_" + filename_t + "_" + filename_t2 + ".pickle"
        fitname = "cro_" + filename_t + "_" + filename_t2 + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['CRO'] = roc_auc

        ###############################

        # CEFTAZIDIME

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['CAZ']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'CAZ', y)

        # set Y variable
        urines5['Y'] = urines5['CAZ']
        urines5 = urines5.drop('CAZ', axis=1)
        features = urines5.drop('Y', axis=1)

        urines6 = pd.read_csv(filename_t2)
        urines6['standard_age'] = urines6['standard_age'].map(str)
        y = urines6['CAZ']
        urines6 = pd.get_dummies(urines6.drop(drops, axis=1))
        urines6.insert(0, 'CAZ', y)

        # set Y variable
        urines6['Y'] = urines6['CAZ']
        urines6 = urines6.drop('CAZ', axis=1)
        features2 = urines6.drop('Y', axis=1)

        missing_columns = set(features.columns) - set(features2.columns)
        for col in missing_columns:
            features2[col] = False

        features2 = features2[features.columns]

        # C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                          target=urines5["Y"],
                          ht_features=features,
                          test=0.2,
                          random=i,
                          folds=6,
                          model_type="ovr",
                          cw='balanced',
                          scorer='roc_auc_ovr',
                          target_ab='Ceftazidime')

        vari_overall['CAZ'] = vari_list
        c_values['CAZ'] = c_value

        # LR final run
        LR_multi_final_t(target=urines5["Y"],
                       final_features=features,
                       test=0.2,
                       random=i,
                       C_val=c_value,
                       cw='balanced',
                       model_type="ovr",
                       thr=0.5,
                       ab_of_interest="Ceftazidime",
                       av="micro",
                       target2=urines6["Y"],
                       test_datf=features2)

        class_reps['CAZ'] = class_report

        # save fit, variables, and hyperparameters
        filename = "caz_" + filename_t + "_" + filename_t2 + ".pickle"
        fitname = "caz_" + filename_t + "_" + filename_t2 + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['CAZ'] = roc_auc

        ###############################

        # CEFEPIME

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['FEP']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'FEP', y)

        # set Y variable
        urines5['Y'] = urines5['FEP']
        urines5 = urines5.drop('FEP', axis=1)
        features = urines5.drop('Y', axis=1)

        urines6 = pd.read_csv(filename_t2)
        urines6['standard_age'] = urines6['standard_age'].map(str)
        y = urines6['FEP']
        urines6 = pd.get_dummies(urines6.drop(drops, axis=1))
        urines6.insert(0, 'FEP', y)

        # set Y variable
        urines6['Y'] = urines6['FEP']
        urines6 = urines6.drop('FEP', axis=1)
        features2 = urines6.drop('Y', axis=1)

        missing_columns = set(features.columns) - set(features2.columns)
        for col in missing_columns:
            features2[col] = False

        features2 = features2[features.columns]

        # C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                          target=urines5["Y"],
                          ht_features=features,
                          test=0.2,
                          random=i,
                          folds=6,
                          model_type="ovr",
                          cw='balanced',
                          scorer='roc_auc_ovr',
                          target_ab='Cefepime')

        vari_overall['FEP'] = vari_list
        c_values['FEP'] = c_value

        # LR final run
        LR_multi_final_t(target=urines5["Y"],
                       final_features=features,
                       test=0.2,
                       random=i,
                       C_val=c_value,
                       cw='balanced',
                       model_type="ovr",
                       thr=0.5,
                       ab_of_interest="Cefepime",
                       av="micro",
                       target2=urines6["Y"],
                       test_datf=features2)

        class_reps['FEP'] = class_report

        # save fit, variables, and hyperparameters
        filename = "fep_" + filename_t + "_" + filename_t2 + ".pickle"
        fitname = "fep_" + filename_t + "_" + filename_t2 + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['FEP'] = roc_auc

        ###############################

        # MEROPENEM

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['MEM']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'MEM', y)

        # set Y variable
        urines5['Y'] = urines5['MEM']
        urines5 = urines5.drop('MEM', axis=1)
        features = urines5.drop('Y', axis=1)

        urines6 = pd.read_csv(filename_t2)
        urines6['standard_age'] = urines6['standard_age'].map(str)
        y = urines6['MEM']
        urines6 = pd.get_dummies(urines6.drop(drops, axis=1))
        urines6.insert(0, 'MEM', y)

        # set Y variable
        urines6['Y'] = urines6['MEM']
        urines6 = urines6.drop('MEM', axis=1)
        features2 = urines6.drop('Y', axis=1)

        missing_columns = set(features.columns) - set(features2.columns)
        for col in missing_columns:
            features2[col] = False

        features2 = features2[features.columns]

        # C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                          target=urines5["Y"],
                          ht_features=features,
                          test=0.2,
                          random=i,
                          folds=6,
                          model_type="ovr",
                          cw='balanced',
                          scorer='roc_auc_ovr',
                          target_ab='Meropenem')

        vari_overall['MEM'] = vari_list
        c_values['MEM'] = c_value

        # LR final run
        LR_multi_final_t(target=urines5["Y"],
                       final_features=features,
                       test=0.2,
                       random=i,
                       C_val=c_value,
                       cw='balanced',
                       model_type="ovr",
                       thr=0.5,
                       ab_of_interest="Meropenem",
                       av="micro",
                       target2=urines6["Y"],
                       test_datf=features2)

        class_reps['MEM'] = class_report

        # save fit, variables, and hyperparameters
        filename = "mem_" + filename_t + "_" + filename_t2 + ".pickle"
        fitname = "mem_" + filename_t + "_" + filename_t2 + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['MEM'] = roc_auc

        ###############################

        # CIPROFLOXACIN

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['CIP']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'CIP', y)

        # set Y variable
        urines5['Y'] = urines5['CIP']
        urines5 = urines5.drop('CIP', axis=1)
        features = urines5.drop('Y', axis=1)

        urines6 = pd.read_csv(filename_t2)
        urines6['standard_age'] = urines6['standard_age'].map(str)
        y = urines6['CIP']
        urines6 = pd.get_dummies(urines6.drop(drops, axis=1))
        urines6.insert(0, 'CIP', y)

        # set Y variable
        urines6['Y'] = urines6['CIP']
        urines6 = urines6.drop('CIP', axis=1)
        features2 = urines6.drop('Y', axis=1)

        missing_columns = set(features.columns) - set(features2.columns)
        for col in missing_columns:
            features2[col] = False

        features2 = features2[features.columns]

        # C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                          target=urines5["Y"],
                          ht_features=features,
                          test=0.2,
                          random=i,
                          folds=6,
                          model_type="ovr",
                          cw='balanced',
                          scorer='roc_auc_ovr',
                          target_ab='Ciprofloxacin')

        vari_overall['CIP'] = vari_list
        c_values['CIP'] = c_value

        # LR final run
        LR_multi_final_t(target=urines5["Y"],
                       final_features=features,
                       test=0.2,
                       random=i,
                       C_val=c_value,
                       cw='balanced',
                       model_type="ovr",
                       thr=0.5,
                       ab_of_interest="Ciprofloxacin",
                       av="micro",
                       target2=urines6["Y"],
                       test_datf=features2)

        class_reps['CIP'] = class_report

        # save fit, variables, and hyperparameters
        filename = "cip_" + filename_t + "_" + filename_t2 + ".pickle"
        fitname = "cip_" + filename_t + "_" + filename_t2 + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['CIP'] = roc_auc

        ###############################

        # GENTAMICIN

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['GEN']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'GEN', y)

        # set Y variable
        urines5['Y'] = urines5['GEN']
        urines5 = urines5.drop('GEN', axis=1)
        features = urines5.drop('Y', axis=1)

        urines6 = pd.read_csv(filename_t2)
        urines6['standard_age'] = urines6['standard_age'].map(str)
        y = urines6['GEN']
        urines6 = pd.get_dummies(urines6.drop(drops, axis=1))
        urines6.insert(0, 'GEN', y)

        # set Y variable
        urines6['Y'] = urines6['GEN']
        urines6 = urines6.drop('GEN', axis=1)
        features2 = urines6.drop('Y', axis=1)

        missing_columns = set(features.columns) - set(features2.columns)
        for col in missing_columns:
            features2[col] = False

        features2 = features2[features.columns]

        # C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                          target=urines5["Y"],
                          ht_features=features,
                          test=0.2,
                          random=i,
                          folds=6,
                          model_type="ovr",
                          cw='balanced',
                          scorer='roc_auc_ovr',
                          target_ab='Gentamicin')

        vari_overall['GEN'] = vari_list
        c_values['GEN'] = c_value

        # LR final run
        LR_multi_final_t(target=urines5["Y"],
                       final_features=features,
                       test=0.2,
                       random=i,
                       C_val=c_value,
                       cw='balanced',
                       model_type="ovr",
                       thr=0.5,
                       ab_of_interest="Gentamicin",
                       av="micro",
                       target2=urines6["Y"],
                       test_datf=features2)

        class_reps['GEN'] = class_report

        # save fit, variables, and hyperparameters
        filename = "gen_" + filename_t + "_" + filename_t2 + ".pickle"
        fitname = "gen_" + filename_t + "_" + filename_t2 + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['GEN'] = roc_auc

        ###############################

        # TRIMETHOPRIM-SULFAMETHOXAZOLE

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['SXT']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'SXT', y)

        # set Y variable
        urines5['Y'] = urines5['SXT']
        urines5 = urines5.drop('SXT', axis=1)
        features = urines5.drop('Y', axis=1)

        urines6 = pd.read_csv(filename_t2)
        urines6['standard_age'] = urines6['standard_age'].map(str)
        y = urines6['SXT']
        urines6 = pd.get_dummies(urines6.drop(drops, axis=1))
        urines6.insert(0, 'SXT', y)

        # set Y variable
        urines6['Y'] = urines6['SXT']
        urines6 = urines6.drop('SXT', axis=1)
        features2 = urines6.drop('Y', axis=1)

        missing_columns = set(features.columns) - set(features2.columns)
        for col in missing_columns:
            features2[col] = False

        features2 = features2[features.columns]

        # C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                          target=urines5["Y"],
                          ht_features=features,
                          test=0.2,
                          random=i,
                          folds=6,
                          model_type="ovr",
                          cw='balanced',
                          scorer='roc_auc_ovr',
                          target_ab='Trimethoprim-sulfamethoxazole',
                          targ_result='R')

        vari_overall['SXT'] = vari_list
        c_values['SXT'] = c_value

        # LR final run
        LR_multi_final_t(target=urines5["Y"],
                       final_features=features,
                       test=0.2,
                       random=i,
                       C_val=c_value,
                       cw='balanced',
                       model_type="ovr",
                       thr=0.5,
                       ab_of_interest="Trimethoprim-sulfamethoxazole",
                       av="micro",
                       target2=urines6["Y"],
                       test_datf=features2)

        class_reps['SXT'] = class_report

        # save fit, variables, and hyperparameters
        filename = "sxt_" + filename_t + "_" + filename_t2 + ".pickle"
        fitname = "sxt_" + filename_t + "_" + filename_t2 + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['SXT'] = roc_auc

        ###############################

        # NITROFURANTOIN

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['NIT']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'NIT', y)

        # set Y variable
        urines5['Y'] = urines5['NIT']
        urines5 = urines5.drop('NIT', axis=1)
        features = urines5.drop('Y', axis=1)

        urines6 = pd.read_csv(filename_t2)
        urines6['standard_age'] = urines6['standard_age'].map(str)
        y = urines6['NIT']
        urines6 = pd.get_dummies(urines6.drop(drops, axis=1))
        urines6.insert(0, 'NIT', y)

        # set Y variable
        urines6['Y'] = urines6['NIT']
        urines6 = urines6.drop('NIT', axis=1)
        features2 = urines6.drop('Y', axis=1)

        missing_columns = set(features.columns) - set(features2.columns)
        for col in missing_columns:
            features2[col] = False

        features2 = features2[features.columns]

        # C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                          target=urines5["Y"],
                          ht_features=features,
                          test=0.2,
                          random=i,
                          folds=6,
                          model_type="ovr",
                          cw='balanced',
                          scorer='roc_auc_ovr',
                          target_ab='Nitrofurantoin')

        vari_overall['NIT'] = vari_list
        c_values['NIT'] = c_value

        # LR final run
        LR_multi_final_t(target=urines5["Y"],
                       final_features=features,
                       test=0.2,
                       random=i,
                       C_val=c_value,
                       cw='balanced',
                       model_type="ovr",
                       thr=0.5,
                       ab_of_interest="Nitrofurantoin",
                       av="micro",
                       target2=urines6["Y"],
                       test_datf=features2)

        class_reps['NIT'] = class_report

        # save fit, variables, and hyperparameters
        filename = "nit_" + filename_t + "_" + filename_t2 + ".pickle"
        fitname = "nit_" + filename_t + "_" + filename_t2 + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['NIT'] = roc_auc

        ###############################

        # VANCOMYCIN

        ###############################

        urines5 = pd.read_csv(filename_t)
        urines5['standard_age'] = urines5['standard_age'].map(str)
        y = urines5['VAN']
        urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
        urines5.insert(0, 'VAN', y)

        # set Y variable
        urines5['Y'] = urines5['VAN']
        urines5 = urines5.drop('VAN', axis=1)
        features = urines5.drop('Y', axis=1)

        urines6 = pd.read_csv(filename_t2)
        urines6['standard_age'] = urines6['standard_age'].map(str)
        y = urines6['VAN']
        urines6 = pd.get_dummies(urines6.drop(drops, axis=1))
        urines6.insert(0, 'VAN', y)

        # set Y variable
        urines6['Y'] = urines6['VAN']
        urines6 = urines6.drop('VAN', axis=1)
        features2 = urines6.drop('Y', axis=1)

        missing_columns = set(features.columns) - set(features2.columns)
        for col in missing_columns:
            features2[col] = False

        features2 = features2[features.columns]

        # C hyperparameter tuning and feature selection
        lr_hyp_tune_feats_t(df=urines5,
                          target=urines5["Y"],
                          ht_features=features,
                          test=0.2,
                          random=i,
                          folds=6,
                          model_type="ovr",
                          cw='balanced',
                          scorer='roc_auc_ovr',
                          target_ab='Vancomycin',
                          targ_result='R')

        vari_overall['VAN'] = vari_list
        c_values['VAN'] = c_value

        # LR final run
        LR_multi_final_t(target=urines5["Y"],
                       final_features=features,
                       test=0.2,
                       random=i,
                       C_val=c_value,
                       cw='balanced',
                       model_type="ovr",
                       thr=0.5,
                       ab_of_interest="Vancomycin",
                       av="micro",
                       target2=urines6["Y"],
                       test_datf=features2)

        class_reps['VAN'] = class_report

        # save fit, variables, and hyperparameters
        filename = "van_" + filename_t + "_" + filename_t2 + ".pickle"
        fitname = "van_" + filename_t + "_" + filename_t2 + ".pickle"
        with open(filename, 'wb') as f:
            pickle.dump(vari_list, f)
        with open(fitname, 'wb') as f:
            pickle.dump(model_dict, f)

        rocs['VAN'] = roc_auc

        aucrocs[i] = rocs

across_time_per("urines_5t_1.csv","urines_5t_2.csv")
aucrocs_5t_1_2 = aucrocs
aucrocs_5t_1_2['df1'] = 1
aucrocs_5t_1_2['df2'] = 2
aucrocs_5t_1_2 = pd.DataFrame(aucrocs_5t_1_2)

across_time_per("urines_5t_1.csv","urines_5t_3.csv")
aucrocs_5t_1_3 = aucrocs
aucrocs_5t_1_3['df1'] = 1
aucrocs_5t_1_3['df2'] = 3
aucrocs_5t_1_3 = pd.DataFrame(aucrocs_5t_1_3)

across_time_per("urines_5t_1.csv","urines_5t_4.csv")
aucrocs_5t_1_4 = aucrocs
aucrocs_5t_1_4['df1'] = 1
aucrocs_5t_1_4['df2'] = 4
aucrocs_5t_1_4 = pd.DataFrame(aucrocs_5t_1_4)

across_time_per("urines_5t_2.csv","urines_5t_1.csv")
aucrocs_5t_2_1 = aucrocs
aucrocs_5t_2_1['df1'] = 2
aucrocs_5t_2_1['df2'] = 1
aucrocs_5t_2_1 = pd.DataFrame(aucrocs_5t_2_1)

across_time_per("urines_5t_2.csv","urines_5t_3.csv")
aucrocs_5t_2_3 = aucrocs
aucrocs_5t_2_3['df1'] = 2
aucrocs_5t_2_3['df2'] = 3
aucrocs_5t_2_3 = pd.DataFrame(aucrocs_5t_2_3)

across_time_per("urines_5t_2.csv","urines_5t_4.csv")
aucrocs_5t_2_4 = aucrocs
aucrocs_5t_2_4['df1'] = 2
aucrocs_5t_2_4['df2'] = 4
aucrocs_5t_2_4 = pd.DataFrame(aucrocs_5t_2_4)

across_time_per("urines_5t_3.csv","urines_5t_1.csv")
aucrocs_5t_3_1 = aucrocs
aucrocs_5t_3_1['df1'] = 3
aucrocs_5t_3_1['df2'] = 1
aucrocs_5t_3_1 = pd.DataFrame(aucrocs_5t_3_1)

across_time_per("urines_5t_3.csv","urines_5t_2.csv")
aucrocs_5t_3_2 = aucrocs
aucrocs_5t_3_2['df1'] = 3
aucrocs_5t_3_2['df2'] = 2
aucrocs_5t_3_2 = pd.DataFrame(aucrocs_5t_3_2)

across_time_per("urines_5t_3.csv","urines_5t_4.csv")
aucrocs_5t_3_4 = aucrocs
aucrocs_5t_3_4['df1'] = 3
aucrocs_5t_3_4['df2'] = 4
aucrocs_5t_3_4 = pd.DataFrame(aucrocs_5t_3_4)

across_time_per("urines_5t_4.csv","urines_5t_1.csv")
aucrocs_5t_4_1 = aucrocs
aucrocs_5t_4_1['df1'] = 4
aucrocs_5t_4_1['df2'] = 1
aucrocs_5t_4_1 = pd.DataFrame(aucrocs_5t_4_1)

across_time_per("urines_5t_4.csv","urines_5t_2.csv")
aucrocs_5t_4_2 = aucrocs
aucrocs_5t_4_2['df1'] = 4
aucrocs_5t_4_2['df2'] = 2
aucrocs_5t_4_2 = pd.DataFrame(aucrocs_5t_4_2)

across_time_per("urines_5t_4.csv","urines_5t_3.csv")
aucrocs_5t_4_3 = aucrocs
aucrocs_5t_4_3['df1'] = 4
aucrocs_5t_4_3['df2'] = 3
aucrocs_5t_4_3 = pd.DataFrame(aucrocs_5t_4_3)

aucrocs_5t = pd.concat([aucrocs_5t_1_1, aucrocs_5t_1_2,
                     aucrocs_5t_1_3,aucrocs_5t_1_4,
                     aucrocs_5t_2_1,aucrocs_5t_2_2,
                     aucrocs_5t_2_3,aucrocs_5t_2_4,
                     aucrocs_5t_3_1,aucrocs_5t_3_2,
                     aucrocs_5t_3_3,aucrocs_5t_3_4,
                     aucrocs_5t_4_1,aucrocs_5t_4_2,
                     aucrocs_5t_4_3,aucrocs_5t_4_4],ignore_index=True)

aucrocs_5t.to_csv('aucrocs_5t.csv',index=False)
