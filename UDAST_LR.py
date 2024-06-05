#load in, check and clean

import pandas as pd

drops = ["AMP","SAM","TZP","CZO","CRO","CAZ","FEP","MEM","CIP","GEN","SXT","NIT","VAN" ]
ref_df = pd.read_csv("urines5.csv")
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

###############################

#AMPICILLIN

###############################

urines5 = pd.read_csv("urines5.csv")
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

##protected variable analysis

amp_fair_classrepsT = {}
amp_fair_classrepsF = {}

for protvar in prot_vars:

    mod_fairness(protvarv = protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                   final_features=features,
                   test=0.2,
                   random=100,
                   C_val=c_value,
                   cw='balanced',
                   model_type="ovr",
                   thr=0.5,
                   ab_of_interest="Ampicillin",
                   av="micro",
                 Bool=True)

    amp_fair_classrepsT[protvar] = class_report

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Ampicillin",
                 av="micro",
                 Bool=False)

    amp_fair_classrepsF[protvar] = class_report


#Stability prediction
p2 = pd.DataFrame()
p3 = pd.DataFrame()
p4 = pd.DataFrame()
p5 = pd.DataFrame()
p6 = pd.DataFrame()
p7 = pd.DataFrame()
p8 = pd.DataFrame()
p9 = pd.DataFrame()

iter_stab_pred("Ampicillin")

p2_amp = p2
p3_amp = p3
p4_amp = p4
p5_amp = p5
p6_amp = p6
p7_amp = p7
p8_amp = p8
p9_amp = p9

p2_amp.to_csv("p2_amp.csv")
p3_amp.to_csv("p3_amp.csv")
p4_amp.to_csv("p4_amp.csv")
p5_amp.to_csv("p5_amp.csv")
p6_amp.to_csv("p6_amp.csv")
p7_amp.to_csv("p7_amp.csv")
p8_amp.to_csv("p8_amp.csv")
p9_amp.to_csv("p9_amp.csv")


###############################

#AMPICILLIN-SULBACTAM

###############################

urines5 = pd.read_csv("urines5.csv")
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

##protected variable analysis

sam_fair_classrepsT = {}
sam_fair_classrepsF = {}

for protvar in prot_vars:

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Ampicillin-sulbactam",
                 av="micro",
                 Bool=True)

    sam_fair_classrepsT[protvar] = class_report


    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Ampicillin-sulbactam",
                 av="micro",
                 Bool=False)

    sam_fair_classrepsF[protvar] = class_report

#Stability prediction
p2 = pd.DataFrame()
p3 = pd.DataFrame()
p4 = pd.DataFrame()
p5 = pd.DataFrame()
p6 = pd.DataFrame()
p7 = pd.DataFrame()
p8 = pd.DataFrame()
p9 = pd.DataFrame()

iter_stab_pred("Ampicillin-sulbactam")

p2_sam = p2
p3_sam = p3
p4_sam = p4
p5_sam = p5
p6_sam = p6
p7_sam = p7
p8_sam = p8
p9_sam = p9

p2_sam.to_csv("p2_sam.csv")
p3_sam.to_csv("p3_sam.csv")
p4_sam.to_csv("p4_sam.csv")
p5_sam.to_csv("p5_sam.csv")
p6_sam.to_csv("p6_sam.csv")
p7_sam.to_csv("p7_sam.csv")
p8_sam.to_csv("p8_sam.csv")
p9_sam.to_csv("p9_sam.csv")

###############################

#PIPERACILLIN-TAZOBACTAM

###############################

urines5 = pd.read_csv("urines5.csv")
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

##protected variable analysis

tzp_fair_classrepsT = {}
tzp_fair_classrepsF = {}

for protvar in prot_vars:

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Piperacillin-tazobactam",
                 av="micro",
                 Bool=True)

    tzp_fair_classrepsT[protvar] = class_report

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Piperacillin-tazobactam",
                 av="micro",
                 Bool=False)

    tzp_fair_classrepsF[protvar] = class_report

#Stability prediction
p2 = pd.DataFrame()
p3 = pd.DataFrame()
p4 = pd.DataFrame()
p5 = pd.DataFrame()
p6 = pd.DataFrame()
p7 = pd.DataFrame()
p8 = pd.DataFrame()
p9 = pd.DataFrame()

iter_stab_pred("Piperacillin-tazobactam")

p2_tzp = p2
p3_tzp = p3
p4_tzp = p4
p5_tzp = p5
p6_tzp = p6
p7_tzp = p7
p8_tzp = p8
p9_tzp = p9

p2_tzp.to_csv("p2_tzp.csv")
p3_tzp.to_csv("p3_tzp.csv")
p4_tzp.to_csv("p4_tzp.csv")
p5_tzp.to_csv("p5_tzp.csv")
p6_tzp.to_csv("p6_tzp.csv")
p7_tzp.to_csv("p7_tzp.csv")
p8_tzp.to_csv("p8_tzp.csv")
p9_tzp.to_csv("p9_tzp.csv")

###############################

#CEFAZOLIN

###############################

urines5 = pd.read_csv("urines5.csv")
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

##protected variable analysis

czo_fair_classrepsT = {}
czo_fair_classrepsF = {}

for protvar in prot_vars:

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Cefazolin",
                 av="micro",
                 Bool=True)

    czo_fair_classrepsT[protvar] = class_report

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Cefazolin",
                 av="micro",
                 Bool=False)

    czo_fair_classrepsF[protvar] = class_report

#Stability prediction
p2 = pd.DataFrame()
p3 = pd.DataFrame()
p4 = pd.DataFrame()
p5 = pd.DataFrame()
p6 = pd.DataFrame()
p7 = pd.DataFrame()
p8 = pd.DataFrame()
p9 = pd.DataFrame()

iter_stab_pred("Cefazolin")

p2_czo = p2
p3_czo = p3
p4_czo = p4
p5_czo = p5
p6_czo = p6
p7_czo = p7
p8_czo = p8
p9_czo = p9

p2_czo.to_csv("p2_czo.csv")
p3_czo.to_csv("p3_czo.csv")
p4_czo.to_csv("p4_czo.csv")
p5_czo.to_csv("p5_czo.csv")
p6_czo.to_csv("p6_czo.csv")
p7_czo.to_csv("p7_czo.csv")
p8_czo.to_csv("p8_czo.csv")
p9_czo.to_csv("p9_czo.csv")

###############################

#CEFTRIAXONE

###############################

urines5 = pd.read_csv("urines5.csv")
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

##protected variable analysis

cro_fair_classrepsT = {}
cro_fair_classrepsF = {}

for protvar in prot_vars:
    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Ceftriaxone",
                 av="micro",
                 Bool=True)

    cro_fair_classrepsT[protvar] = class_report

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Ceftriaxone",
                 av="micro",
                 Bool=True)

    cro_fair_classrepsF[protvar] = class_report

#Stability prediction
p2 = pd.DataFrame()
p3 = pd.DataFrame()
p4 = pd.DataFrame()
p5 = pd.DataFrame()
p6 = pd.DataFrame()
p7 = pd.DataFrame()
p8 = pd.DataFrame()
p9 = pd.DataFrame()

iter_stab_pred("Ceftriaxone")

p2_cro = p2
p3_cro = p3
p4_cro = p4
p5_cro = p5
p6_cro = p6
p7_cro = p7
p8_cro = p8
p9_cro = p9

p2_cro.to_csv("p2_cro.csv")
p3_cro.to_csv("p3_cro.csv")
p4_cro.to_csv("p4_cro.csv")
p5_cro.to_csv("p5_cro.csv")
p6_cro.to_csv("p6_cro.csv")
p7_cro.to_csv("p7_cro.csv")
p8_cro.to_csv("p8_cro.csv")
p9_cro.to_csv("p9_cro.csv")

###############################

#CEFTAZIDIME

###############################

urines5 = pd.read_csv("urines5.csv")
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

##protected variable analysis

caz_fair_classrepsT = {}
caz_fair_classrepsF = {}

for protvar in prot_vars:
    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Ceftazidime",
                 av="micro",
                 Bool=True)

    caz_fair_classrepsT[protvar] = class_report

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Ceftazidime",
                 av="micro",
                 Bool=False)

    caz_fair_classrepsF[protvar] = class_report

#Stability prediction
p2 = pd.DataFrame()
p3 = pd.DataFrame()
p4 = pd.DataFrame()
p5 = pd.DataFrame()
p6 = pd.DataFrame()
p7 = pd.DataFrame()
p8 = pd.DataFrame()
p9 = pd.DataFrame()

iter_stab_pred("Ceftazidime")

p2_caz = p2
p3_caz = p3
p4_caz = p4
p5_caz = p5
p6_caz = p6
p7_caz = p7
p8_caz = p8
p9_caz = p9

p2_caz.to_csv("p2_caz.csv")
p3_caz.to_csv("p3_caz.csv")
p4_caz.to_csv("p4_caz.csv")
p5_caz.to_csv("p5_caz.csv")
p6_caz.to_csv("p6_caz.csv")
p7_caz.to_csv("p7_caz.csv")
p8_caz.to_csv("p8_caz.csv")
p9_caz.to_csv("p9_caz.csv")

###############################

#CEFEPIME

###############################

urines5 = pd.read_csv("urines5.csv")
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

##protected variable analysis

fep_fair_classrepsT = {}
fep_fair_classrepsF = {}

for protvar in prot_vars:

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Cefepime",
                 av="micro",
                 Bool=True)

    fep_fair_classrepsT[protvar] = class_report

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Cefepime",
                 av="micro",
                 Bool=False)

    fep_fair_classrepsF[protvar] = class_report

#Stability prediction
iter_stab_pred("Cefepime")

p2_fep = p2
p3_fep = p3
p4_fep = p4
p5_fep = p5
p6_fep = p6
p7_fep = p7
p8_fep = p8
p9_fep = p9

p2_fep.to_csv("p2_fep.csv")
p3_fep.to_csv("p3_fep.csv")
p4_fep.to_csv("p4_fep.csv")
p5_fep.to_csv("p5_fep.csv")
p6_fep.to_csv("p6_fep.csv")
p7_fep.to_csv("p7_fep.csv")
p8_fep.to_csv("p8_fep.csv")
p9_fep.to_csv("p9_fep.csv")

###############################

#MEROPENEM

###############################

urines5 = pd.read_csv("urines5.csv")
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

##protected variable analysis

mem_fair_classrepsT = {}
mem_fair_classrepsF = {}

for protvar in prot_vars:
    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Meropenem",
                 av="micro",
                 Bool=True)

    mem_fair_classrepsT[protvar] = class_report

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Meropenem",
                 av="micro",
                 Bool=False)

    mem_fair_classrepsF[protvar] = class_report

#Stability prediction
p2 = pd.DataFrame()
p3 = pd.DataFrame()
p4 = pd.DataFrame()
p5 = pd.DataFrame()
p6 = pd.DataFrame()
p7 = pd.DataFrame()
p8 = pd.DataFrame()
p9 = pd.DataFrame()

iter_stab_pred("Meropenem")

p2_mem = p2
p3_mem = p3
p4_mem = p4
p5_mem = p5
p6_mem = p6
p7_mem = p7
p8_mem = p8
p9_mem = p9

p2_mem.to_csv("p2_mem.csv")
p3_mem.to_csv("p3_mem.csv")
p4_mem.to_csv("p4_mem.csv")
p5_mem.to_csv("p5_mem.csv")
p6_mem.to_csv("p6_mem.csv")
p7_mem.to_csv("p7_mem.csv")
p8_mem.to_csv("p8_mem.csv")
p9_mem.to_csv("p9_mem.csv")

###############################

#CIPROFLOXACIN

###############################

urines5 = pd.read_csv("urines5.csv")
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

##protected variable analysis

cip_fair_classrepsT = {}
cip_fair_classrepsF = {}

for protvar in prot_vars:
    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Ciprofloxacin",
                 av="micro",
                 Bool=False)

    cip_fair_classrepsT[protvar] = class_report

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Ciprofloxacin",
                 av="micro",
                 Bool=False)

    cip_fair_classrepsF[protvar] = class_report

#Stability prediction
p2 = pd.DataFrame()
p3 = pd.DataFrame()
p4 = pd.DataFrame()
p5 = pd.DataFrame()
p6 = pd.DataFrame()
p7 = pd.DataFrame()
p8 = pd.DataFrame()
p9 = pd.DataFrame()

iter_stab_pred("Ciprofloxacin")

p2_cip = p2
p3_cip = p3
p4_cip = p4
p5_cip = p5
p6_cip = p6
p7_cip = p7
p8_cip = p8
p9_cip = p9

p2_cip.to_csv("p2_cip.csv")
p3_cip.to_csv("p3_cip.csv")
p4_cip.to_csv("p4_cip.csv")
p5_cip.to_csv("p5_cip.csv")
p6_cip.to_csv("p6_cip.csv")
p7_cip.to_csv("p7_cip.csv")
p8_cip.to_csv("p8_cip.csv")
p9_cip.to_csv("p9_cip.csv")

###############################

#GENTAMICIN

###############################

urines5 = pd.read_csv("urines5.csv")
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

##protected variable analysis

gen_fair_classrepsT = {}
gen_fair_classrepsF = {}

for protvar in prot_vars:

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Gentamicin",
                 av="micro",
                 Bool=True)

    gen_fair_classrepsT[protvar] = class_report

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Gentamicin",
                 av="micro",
                 Bool=False)

    gen_fair_classrepsF[protvar] = class_report

#Stability prediction
p2 = pd.DataFrame()
p3 = pd.DataFrame()
p4 = pd.DataFrame()
p5 = pd.DataFrame()
p6 = pd.DataFrame()
p7 = pd.DataFrame()
p8 = pd.DataFrame()
p9 = pd.DataFrame()

iter_stab_pred("Gentamicin")

p2_gen = p2
p3_gen = p3
p4_gen = p4
p5_gen = p5
p6_gen = p6
p7_gen = p7
p8_gen = p8
p9_gen = p9

p2_gen.to_csv("p2_gen.csv")
p3_gen.to_csv("p3_gen.csv")
p4_gen.to_csv("p4_gen.csv")
p5_gen.to_csv("p5_gen.csv")
p6_gen.to_csv("p6_gen.csv")
p7_gen.to_csv("p7_gen.csv")
p8_gen.to_csv("p8_gen.csv")
p9_gen.to_csv("p9_gen.csv")

###############################

#TRIMETHOPRIM-SULFAMETHOXAZOLE

###############################

urines5 = pd.read_csv("urines5.csv")
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

##protected variable analysis

sxt_fair_classrepsT = {}
sxt_fair_classrepsF = {}

for protvar in prot_vars:

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Trimethoprim-sulfamethoxazole",
                 av="micro",
                 Bool=True)

    sxt_fair_classrepsT[protvar] = class_report

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Trimethoprim-sulfamethoxazole",
                 av="micro",
                 Bool=False)

    sxt_fair_classrepsF[protvar] = class_report

#Stability prediction
p2 = pd.DataFrame()
p3 = pd.DataFrame()
p4 = pd.DataFrame()
p5 = pd.DataFrame()
p6 = pd.DataFrame()
p7 = pd.DataFrame()
p8 = pd.DataFrame()
p9 = pd.DataFrame()

iter_stab_pred("Trimethoprim-sulfamethoxazole")

p2_sxt = p2
p3_sxt = p3
p4_sxt = p4
p5_sxt = p5
p6_sxt = p6
p7_sxt = p7
p8_sxt = p8
p9_sxt = p9

p2_sxt.to_csv("p2_sxt.csv")
p3_sxt.to_csv("p3_sxt.csv")
p4_sxt.to_csv("p4_sxt.csv")
p5_sxt.to_csv("p5_sxt.csv")
p6_sxt.to_csv("p6_sxt.csv")
p7_sxt.to_csv("p7_sxt.csv")
p8_sxt.to_csv("p8_sxt.csv")
p9_sxt.to_csv("p9_sxt.csv")

###############################

#NITROFURANTOIN

###############################

urines5 = pd.read_csv("urines5.csv")
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

##protected variable analysis

nit_fair_classrepsT = {}
nit_fair_classrepsF = {}

for protvar in prot_vars:

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Nitrofurantoin",
                 av="micro",
                 Bool=True)

    nit_fair_classrepsT[protvar] = class_report

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Nitrofurantoin",
                 av="micro",
                 Bool=False)

    nit_fair_classrepsF[protvar] = class_report

#Stability prediction
p2 = pd.DataFrame()
p3 = pd.DataFrame()
p4 = pd.DataFrame()
p5 = pd.DataFrame()
p6 = pd.DataFrame()
p7 = pd.DataFrame()
p8 = pd.DataFrame()
p9 = pd.DataFrame()

iter_stab_pred("Nitrofurantoin")

p2_nit = p2
p3_nit = p3
p4_nit = p4
p5_nit = p5
p6_nit = p6
p7_nit = p7
p8_nit = p8
p9_nit = p9

p2_nit.to_csv("p2_nit.csv")
p3_nit.to_csv("p3_nit.csv")
p4_nit.to_csv("p4_nit.csv")
p5_nit.to_csv("p5_nit.csv")
p6_nit.to_csv("p6_nit.csv")
p7_nit.to_csv("p7_nit.csv")
p8_nit.to_csv("p8_nit.csv")
p9_nit.to_csv("p9_nit.csv")

###############################

#VANCOMYCIN

###############################

urines5 = pd.read_csv("urines5.csv")
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
##protected variable analysis

van_fair_classrepsT = {}
van_fair_classrepsF = {}

for protvar in prot_vars:

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Vancomycin",
                 av="micro",
                 Bool=True)

    van_fair_classrepsT[protvar] = class_report

    mod_fairness(protvarv=protvar,
                 chosen_model=fairness_model,
                 target=urines5["Y"],
                 final_features=features,
                 test=0.2,
                 random=100,
                 C_val=c_value,
                 cw='balanced',
                 model_type="ovr",
                 thr=0.5,
                 ab_of_interest="Vancomycin",
                 av="micro",
                 Bool=False)

    van_fair_classrepsF[protvar] = class_report

#Stability prediction
p2 = pd.DataFrame()
p3 = pd.DataFrame()
p4 = pd.DataFrame()
p5 = pd.DataFrame()
p6 = pd.DataFrame()
p7 = pd.DataFrame()
p8 = pd.DataFrame()
p9 = pd.DataFrame()

iter_stab_pred("Vancomycin")

p2_van = p2
p3_van = p3
p4_van = p4
p5_van = p5
p6_van = p6
p7_van = p7
p8_van = p8
p9_van = p9

p2_van.to_csv("p2_van.csv")
p3_van.to_csv("p3_van.csv")
p4_van.to_csv("p4_van.csv")
p5_van.to_csv("p5_van.csv")
p6_van.to_csv("p6_van.csv")
p7_van.to_csv("p7_van.csv")
p8_van.to_csv("p8_van.csv")
p9_van.to_csv("p9_van.csv")

with open('classreps.pickle','wb') as f:
    pickle.dump(class_reps,f)
with open('ampfairF.pickle','wb') as f:
    pickle.dump(amp_fair_classrepsF,f)
with open('ampfairT.pickle', 'wb') as f:
    pickle.dump(amp_fair_classrepsT, f)
with open('samfairF.pickle','wb') as f:
    pickle.dump(amp_fair_classrepsF,f)
with open('samfairT.pickle', 'wb') as f:
    pickle.dump(amp_fair_classrepsT, f)
with open('tzpfairF.pickle','wb') as f:
    pickle.dump(amp_fair_classrepsF,f)
with open('tzpfairT.pickle', 'wb') as f:
    pickle.dump(amp_fair_classrepsT, f)
with open('czofairF.pickle','wb') as f:
    pickle.dump(amp_fair_classrepsF,f)
with open('czofairT.pickle', 'wb') as f:
    pickle.dump(amp_fair_classrepsT, f)
with open('crofairF.pickle','wb') as f:
    pickle.dump(amp_fair_classrepsF,f)
with open('crofairT.pickle', 'wb') as f:
    pickle.dump(amp_fair_classrepsT, f)
with open('cazfairF.pickle','wb') as f:
    pickle.dump(amp_fair_classrepsF,f)
with open('cazfairT.pickle', 'wb') as f:
    pickle.dump(amp_fair_classrepsT, f)
with open('fepfairF.pickle','wb') as f:
    pickle.dump(amp_fair_classrepsF,f)
with open('fepfairT.pickle', 'wb') as f:
    pickle.dump(amp_fair_classrepsT, f)
with open('memfairF.pickle','wb') as f:
    pickle.dump(amp_fair_classrepsF,f)
with open('memfairT.pickle', 'wb') as f:
    pickle.dump(amp_fair_classrepsT, f)
with open('cipfairF.pickle','wb') as f:
    pickle.dump(amp_fair_classrepsF,f)
with open('cipfairT.pickle', 'wb') as f:
    pickle.dump(amp_fair_classrepsT, f)
with open('genfairF.pickle','wb') as f:
    pickle.dump(amp_fair_classrepsF,f)
with open('genfairT.pickle', 'wb') as f:
    pickle.dump(amp_fair_classrepsT, f)
with open('sxtfairF.pickle','wb') as f:
    pickle.dump(amp_fair_classrepsF,f)
with open('sxtfairT.pickle', 'wb') as f:
    pickle.dump(amp_fair_classrepsT, f)
with open('nitfairF.pickle','wb') as f:
    pickle.dump(amp_fair_classrepsF,f)
with open('nitfairT.pickle', 'wb') as f:
    pickle.dump(amp_fair_classrepsT, f)
with open('vanfairF.pickle','wb') as f:
    pickle.dump(amp_fair_classrepsF,f)
with open('vanfairT.pickle', 'wb') as f:
    pickle.dump(amp_fair_classrepsT, f)


#retrieve other performance metrics
ab_list = ['AMP','SAM','TZP','CZO','CRO','CAZ','FEP','MEM','CIP','GEN','SXT','NIT']

def result_printer(metric):
    for ab in ab_list:
        ab_r_metric = class_reps[ab][metric]

        print(round(ab_r_metric,2))

#insert metric name as argument
result_printer('accuracy')



#retrieve fairness results

def fairness_printer(df,metric_num):
    for prot in prot_vars:
        print(round(df[prot][metric_num],2))

#insert fairness dataframe as argument
fairness_printer(nit_fair_classrepsT,4)


#coefficients

#insert model files as arguments
with open("nitfile.pickle", 'rb') as f:
    vari_list = pickle.load(f)
with open("nitfits.pickle", 'rb') as f:
    model_dict = pickle.load(f)

print(len(vari_list))
for var in vari_list:
    print(var)

print(model_dict.intercept_)

for coef in model_dict.coef_[0]:
    print(coef)