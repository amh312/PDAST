#AST AVAILABILITY SENSITIVITY ANALYSIS MODEL DEVELOPMENT

##Load in, reference lists, and preprocessing

###Initialising reference objects for selecting relevant columns/indexes in functions
drops1 = ["AMP","SAM","TZP","CZO","CRO","CAZ","FEP","MEM","CIP","GEN","SXT","NIT"]
lc_abs = ["amp","sam","tzp","czo","cro","caz","fep","mem","cip","gen","sxt","nit"]
ab_full = ['Ampicillin','Ampicillin-sulbactam','Piperacillin-tazobactam','Cefazolin','Ceftriaxone',
           'Ceftazidime','Cefepime','Meropenem','Ciprofloxacin','Gentamicin',
           'Trimethoprim-sulfamethoxazole','Nitrofurantoin']
to_drop = ["subject_id","micro_specimen_id"]
metrics_fairness = ['Accuracy','TPR','FPR','FNR','PPP']

###Preprocessing
ref_df = pd.read_csv("urines5b.csv")
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
aucrocs = {}

##Main model development, fairness analysis, and stability analysis

for antimicrobial in list(range(0, len(drops1), 1)):

    ###Training, hyperparameter tuning, and validation

    drops = drops1[antimicrobial]

    train_test_validate_main(csv="urines5_"+lc_abs[antimicrobial]+".csv",
                             abx=drops1[antimicrobial],
                             antibiotic=ab_full[antimicrobial],
                             ab=lc_abs[antimicrobial],
                             type='ast_')
    aucrocs[drops1[antimicrobial]] = roc_auc

##Compiling and exporting results

###Main analysis performance metrics
metrics_df = result_compiler(class_reps,"ast_main_analysis_metrics.csv")
aucrocs_df = pd.DataFrame(list(aucrocs.items()), columns=['Antimicrobial', 'AUC_ROC'])
aucrocs_df['Information'] = 'AST'
aucrocs_df.to_csv("ast_aucrocs.csv", index=False)
