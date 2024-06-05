#load in, check and clean

import pandas as pd

drops = ["AMP","SAM","TZP","CZO","CRO","CAZ","FEP","MEM","CIP","GEN","SXT","NIT","VAN" ]
ref_df = pd.read_csv("urines5c.csv")
ref_df['standard_age'] = ref_df['standard_age'].map(str)
dummy_df = pd.get_dummies(ref_df.drop(drops, axis=1))
vari_overall = {}
c_values = {}
to_drop = ["subject_id","micro_specimen_id"]
class_reps = {}


###############################

#AMPICILLIN

###############################

urines5 = pd.read_csv("urines5c.csv")
urines5['standard_age'] = urines5['standard_age'].map(str)
y = urines5['AMP']
urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
urines5.insert(0, 'AMP', y)


#set Y variable
urines5['Y'] = urines5['AMP']
urines5 = urines5.drop('AMP',axis=1)
features = urines5.drop('Y',axis=1)


#C hyperparameter tuning and feature selection
lr_hyp_tune_feats(df=urines5,
            target=urines5["Y"],
            ht_features=features,
            test=0.2,
            random=100,
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
         random=100,
         C_val=c_value,
         cw='balanced',
         model_type="ovr",
         thr=0.5,
         ab_of_interest="Ampicillin",
         av="micro")

class_reps['AMP'] = class_report

#save fit, variables, and hyperparameters
with open('ampfile.pickle','wb') as f:
    pickle.dump(vari_list,f)
with open('ampfits.pickle','wb') as f:
    pickle.dump(model_dict,f)


###############################

#AMPICILLIN-SULBACTAM

###############################

urines5 = pd.read_csv("urines5c.csv")
urines5['standard_age'] = urines5['standard_age'].map(str)
y = urines5['SAM']
urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
urines5.insert(0, 'SAM', y)

#set Y variable
urines5['Y'] = urines5['SAM']
urines5 = urines5.drop('SAM',axis=1)
features = urines5.drop('Y',axis=1)

#C hyperparameter tuning and feature selection
lr_hyp_tune_feats(df=urines5,
            target=urines5["Y"],
            ht_features=features,
            test=0.2,
            random=100,
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
         random=100,
         C_val=c_value,
         cw='balanced',
         model_type="ovr",
         thr=0.5,
         ab_of_interest="Ampicillin-sulbactam",
         av="micro")

class_reps['SAM'] = class_report

#save fit, variables, and hyperparameters
with open('samfile.pickle','wb') as f:
    pickle.dump(vari_list,f)
with open('samfits.pickle','wb') as f:
    pickle.dump(model_dict,f)


###############################

#PIPERACILLIN-TAZOBACTAM

###############################

urines5 = pd.read_csv("urines5c.csv")
urines5['standard_age'] = urines5['standard_age'].map(str)
y = urines5['TZP']
urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
urines5.insert(0, 'TZP', y)

#set Y variable
urines5['Y'] = urines5['TZP']
urines5 = urines5.drop('TZP',axis=1)
features = urines5.drop('Y',axis=1)

#C hyperparameter tuning and feature selection
lr_hyp_tune_feats(df=urines5,
            target=urines5["Y"],
            ht_features=features,
            test=0.2,
            random=100,
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
         random=100,
         C_val=c_value,
         cw='balanced',
         model_type="ovr",
         thr=0.5,
         ab_of_interest="Piperacillin-tazobactam",
         av="micro")

class_reps['TZP'] = class_report

#save fit, variables, and hyperparameters
with open('tzpfile.pickle','wb') as f:
    pickle.dump(vari_list,f)
with open('tzpfits.pickle','wb') as f:
    pickle.dump(model_dict,f)


###############################

#CEFAZOLIN

###############################

urines5 = pd.read_csv("urines5c.csv")
urines5['standard_age'] = urines5['standard_age'].map(str)
y = urines5['CZO']
urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
urines5.insert(0, 'CZO', y)

#set Y variable
urines5['Y'] = urines5['CZO']
urines5 = urines5.drop('CZO',axis=1)
features = urines5.drop('Y',axis=1)

#C hyperparameter tuning and feature selection
lr_hyp_tune_feats(df=urines5,
            target=urines5["Y"],
            ht_features=features,
            test=0.2,
            random=100,
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
         random=100,
         C_val=c_value,
         cw='balanced',
         model_type="ovr",
         thr=0.5,
         ab_of_interest="Cefazolin",
         av="micro")

class_reps['CZO'] = class_report

#save fit, variables, and hyperparameters
with open('czofile.pickle','wb') as f:
    pickle.dump(vari_list,f)
with open('czofits.pickle','wb') as f:
    pickle.dump(model_dict,f)



###############################

#CEFTRIAXONE

###############################

urines5 = pd.read_csv("urines5c.csv")
urines5['standard_age'] = urines5['standard_age'].map(str)
y = urines5['CRO']
urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
urines5.insert(0, 'CRO', y)

#set Y variable
urines5['Y'] = urines5['CRO']
urines5 = urines5.drop('CRO',axis=1)
features = urines5.drop('Y',axis=1)

#C hyperparameter tuning and feature selection
lr_hyp_tune_feats(df=urines5,
            target=urines5["Y"],
            ht_features=features,
            test=0.2,
            random=100,
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
         random=100,
         C_val=c_value,
         cw='balanced',
         model_type="ovr",
         thr=0.5,
         ab_of_interest="Ceftriaxone",
         av="micro")

class_reps['CRO'] = class_report

#save fit, variables, and hyperparameters
with open('crofile.pickle','wb') as f:
    pickle.dump(vari_list,f)
with open('crofits.pickle','wb') as f:
    pickle.dump(model_dict,f)


###############################

#CEFTAZIDIME

###############################

urines5 = pd.read_csv("urines5c.csv")
urines5['standard_age'] = urines5['standard_age'].map(str)
y = urines5['CAZ']
urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
urines5.insert(0, 'CAZ', y)

#set Y variable
urines5['Y'] = urines5['CAZ']
urines5 = urines5.drop('CAZ',axis=1)
features = urines5.drop('Y',axis=1)

#C hyperparameter tuning and feature selection
lr_hyp_tune_feats(df=urines5,
            target=urines5["Y"],
            ht_features=features,
            test=0.2,
            random=100,
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
         random=100,
         C_val=c_value,
         cw='balanced',
         model_type="ovr",
         thr=0.5,
         ab_of_interest="Ceftazidime",
         av="micro")


class_reps['CAZ'] = class_report

#save fit, variables, and hyperparameters
with open('cazfile.pickle','wb') as f:
    pickle.dump(vari_list,f)
with open('cazfits.pickle','wb') as f:
    pickle.dump(model_dict,f)


###############################

#CEFEPIME

###############################

urines5 = pd.read_csv("urines5c.csv")
urines5['standard_age'] = urines5['standard_age'].map(str)
y = urines5['FEP']
urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
urines5.insert(0, 'FEP', y)

#set Y variable
urines5['Y'] = urines5['FEP']
urines5 = urines5.drop('FEP',axis=1)
features = urines5.drop('Y',axis=1)

#C hyperparameter tuning and feature selection
lr_hyp_tune_feats(df=urines5,
            target=urines5["Y"],
            ht_features=features,
            test=0.2,
            random=100,
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
         random=100,
         C_val=c_value,
         cw='balanced',
         model_type="ovr",
         thr=0.5,
         ab_of_interest="Cefepime",
         av="micro")

class_reps['FEP'] = class_report

#save fit, variables, and hyperparameters
with open('fepfile.pickle','wb') as f:
    pickle.dump(vari_list,f)
with open('fepfits.pickle','wb') as f:
    pickle.dump(model_dict,f)


###############################

#MEROPENEM

###############################

urines5 = pd.read_csv("urines5c.csv")
urines5['standard_age'] = urines5['standard_age'].map(str)
y = urines5['MEM']
urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
urines5.insert(0, 'MEM', y)

#set Y variable
urines5['Y'] = urines5['MEM']
urines5 = urines5.drop('MEM',axis=1)
features = urines5.drop('Y',axis=1)

#C hyperparameter tuning and feature selection
lr_hyp_tune_feats(df=urines5,
            target=urines5["Y"],
            ht_features=features,
            test=0.2,
            random=100,
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
         random=100,
         C_val=c_value,
         cw='balanced',
         model_type="ovr",
         thr=0.5,
         ab_of_interest="Meropenem",
         av="micro")

class_reps['MEM'] = class_report

#save fit, variables, and hyperparameters
with open('memfile.pickle','wb') as f:
    pickle.dump(vari_list,f)
with open('memfits.pickle','wb') as f:
    pickle.dump(model_dict,f)


###############################

#CIPROFLOXACIN

###############################

urines5 = pd.read_csv("urines5c.csv")
urines5['standard_age'] = urines5['standard_age'].map(str)
y = urines5['CIP']
urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
urines5.insert(0, 'CIP', y)

#set Y variable
urines5['Y'] = urines5['CIP']
urines5 = urines5.drop('CIP',axis=1)
features = urines5.drop('Y',axis=1)

#C hyperparameter tuning and feature selection
lr_hyp_tune_feats(df=urines5,
            target=urines5["Y"],
            ht_features=features,
            test=0.2,
            random=100,
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
         random=100,
         C_val=c_value,
         cw='balanced',
         model_type="ovr",
         thr=0.5,
         ab_of_interest="Ciprofloxacin",
         av="micro")

class_reps['CIP'] = class_report

#save fit, variables, and hyperparameters
with open('cipfile.pickle','wb') as f:
    pickle.dump(vari_list,f)
with open('cipfits.pickle','wb') as f:
    pickle.dump(model_dict,f)


###############################

#GENTAMICIN

###############################

urines5 = pd.read_csv("urines5c.csv")
urines5['standard_age'] = urines5['standard_age'].map(str)
y = urines5['GEN']
urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
urines5.insert(0, 'GEN', y)

#set Y variable
urines5['Y'] = urines5['GEN']
urines5 = urines5.drop('GEN',axis=1)
features = urines5.drop('Y',axis=1)

#C hyperparameter tuning and feature selection
lr_hyp_tune_feats(df=urines5,
            target=urines5["Y"],
            ht_features=features,
            test=0.2,
            random=100,
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
         random=100,
         C_val=c_value,
         cw='balanced',
         model_type="ovr",
         thr=0.5,
         ab_of_interest="Gentamicin",
         av="micro")

class_reps['GEN'] = class_report

#save fit, variables, and hyperparameters
with open('genfile.pickle','wb') as f:
    pickle.dump(vari_list,f)
with open('genfits.pickle','wb') as f:
    pickle.dump(model_dict,f)


###############################

#TRIMETHOPRIM-SULFAMETHOXAZOLE

###############################

urines5 = pd.read_csv("urines5c.csv")
urines5['standard_age'] = urines5['standard_age'].map(str)
y = urines5['SXT']
urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
urines5.insert(0, 'SXT', y)

#set Y variable
urines5['Y'] = urines5['SXT']
urines5 = urines5.drop('SXT',axis=1)
features = urines5.drop('Y',axis=1)

#C hyperparameter tuning and feature selection
lr_hyp_tune_feats(df=urines5,
            target=urines5["Y"],
            ht_features=features,
            test=0.2,
            random=100,
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
         random=100,
         C_val=c_value,
         cw='balanced',
         model_type="ovr",
         thr=0.5,
         ab_of_interest="Trimethoprim-sulfamethoxazole",
         av="micro")

class_reps['SXT'] = class_report

#save fit, variables, and hyperparameters
with open('sxtfile.pickle','wb') as f:
    pickle.dump(vari_list,f)
with open('sxtfits.pickle','wb') as f:
    pickle.dump(model_dict,f)


###############################

#NITROFURANTOIN

###############################

urines5 = pd.read_csv("urines5c.csv")
urines5['standard_age'] = urines5['standard_age'].map(str)
y = urines5['NIT']
urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
urines5.insert(0, 'NIT', y)

#set Y variable
urines5['Y'] = urines5['NIT']
urines5 = urines5.drop('NIT',axis=1)
features = urines5.drop('Y',axis=1)

#C hyperparameter tuning and feature selection
lr_hyp_tune_feats(df=urines5,
            target=urines5["Y"],
            ht_features=features,
            test=0.2,
            random=100,
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
         random=100,
         C_val=c_value,
         cw='balanced',
         model_type="ovr",
         thr=0.5,
         ab_of_interest="Nitrofurantoin",
         av="micro")

class_reps['NIT'] = class_report

#save fit, variables, and hyperparameters
with open('nitfile.pickle','wb') as f:
    pickle.dump(vari_list,f)
with open('nitfits.pickle','wb') as f:
    pickle.dump(model_dict,f)


###############################

#VANCOMYCIN

###############################

urines5 = pd.read_csv("urines5c.csv")
urines5['standard_age'] = urines5['standard_age'].map(str)
y = urines5['VAN']
urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
urines5.insert(0, 'VAN', y)

#set Y variable
urines5['Y'] = urines5['VAN']
urines5 = urines5.drop('VAN',axis=1)
features = urines5.drop('Y',axis=1)

#C hyperparameter tuning and feature selection
lr_hyp_tune_feats(df=urines5,
            target=urines5["Y"],
            ht_features=features,
            test=0.2,
            random=100,
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
         random=100,
         C_val=c_value,
         cw='balanced',
         model_type="ovr",
         thr=0.5,
         ab_of_interest="Vancomycin",
         av="micro")

class_reps['VAN'] = class_report

#save fit, variables, and hyperparameters
with open('vanfile.pickle','wb') as f:
    pickle.dump(vari_list,f)
with open('vanfits.pickle','wb') as f:
    pickle.dump(model_dict,f)
