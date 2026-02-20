#MAIN MODEL DEVELOPMENT

##Load in, reference lists, and preprocessing

###Initialising reference objects for selecting relevant columns/indexes in functions
drops = ["AMP","SAM","TZP","CZO","CRO","CAZ","FEP","MEM","CIP","GEN","SXT","NIT","VAN" ]
lc_abs = ["amp","sam","tzp","czo","cro","caz","fep","mem","cip","gen","sxt","nit","van" ]
ab_full = ['Ampicillin','Ampicillin-sulbactam','Piperacillin-tazobactam','Cefazolin','Ceftriaxone',
           'Ceftazidime','Cefepime','Meropenem','Ciprofloxacin','Gentamicin',
           'Trimethoprim-sulfamethoxazole','Nitrofurantoin','Vancomycin']
to_drop = ["subject_id","micro_specimen_id"]
metrics_fairness = ['Accuracy','TPR','FPR','FNR','PPP']

###Preprocessing
ref_df = pd.read_csv("#FILEPATH#/urines5b.csv")
ref_df['standard_age'] = ref_df['standard_age'].map(str)
dummy_df = pd.get_dummies(ref_df.drop(drops, axis=1))

###Initialising empty dictionaries and lists to populate
vari_overall = {}
c_values = {}
class_reps = {}
prot_vars = []
overallfairT = {}
overallfairF = {}
vari_set = {}
intercept_set = {}
coef_set = {}
names_set = {}
aucrocs = {}

##Main model development, fairness analysis, and stability analysis

for antimicrobial in list(range(0, len(drops), 1)):

    ###Training, hyperparameter tuning, and validation
    train_test_validate_main(csv="#FILEPATH#/urines5.csv",
                             abx=drops[antimicrobial],
                             antibiotic=ab_full[antimicrobial],
                             ab=lc_abs[antimicrobial])
    aucrocs[drops[antimicrobial]] = roc_auc

    ###Fairness analysis
    iterate_mod_fairness(df=urines5,
                         feature_df=features,
                         antibiotic=ab_full[antimicrobial],
                         ab=drops[antimicrobial],
                         abx=lc_abs[antimicrobial])

    ###Stability analysis
    iter_stab_pred(df=urines5,
                   feature_df=features,
                   Antimicrobial_agent=ab_full[antimicrobial],
                   abx=lc_abs[antimicrobial])

##Compiling and exporting results

###Main analysis coefficients
urines5.to_csv("#FILEPATH#/urines_coef_df.csv",index=False)

for antimicrobial in list(range(0, len(drops), 1)):
    with open("#FILEPATH#/"+lc_abs[antimicrobial]+"file.pickle", 'rb') as f:
        vari_list = pickle.load(f)
    with open("#FILEPATH#/"+lc_abs[antimicrobial]+"fits.pickle", 'rb') as f:
        model_dict = pickle.load(f)
    df1 = pd.DataFrame({'Parameter': 'Intercept','Value': np.round(model_dict.intercept_,3)})
    df2 = pd.DataFrame({'Parameter': vari_list,'Value': np.round(model_dict.coef_[0],3)})
    df2 = df2.sort_values(by='Value', ascending=False)
    df = pd.concat([df1, df2], ignore_index=True)
    df.to_csv("#FILEPATH#/"+drops[antimicrobial]+'_coef_list.csv',index=False)

###Main analysis performance metrics
metrics_df = result_compiler(class_reps,"#FILEPATH#/main_analysis_metrics.csv")
aucrocs_df = pd.DataFrame(list(aucrocs.items()), columns=['Antimicrobial', 'AUC_ROC'])
aucrocs_df['Information'] = 'Default'
aucrocs_df.to_csv("#FILEPATH#/default_aucrocs.csv", index=False)

###Fairness analysis
fairness_true_metrics = compile_fairness_results(overallfairT,"#FILEPATH#/T_fairness_metrics.csv")
fairness_false_metrics = compile_fairness_results(overallfairF,"#FILEPATH#/F_fairness_metrics.csv")
