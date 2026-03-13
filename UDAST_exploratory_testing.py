
#EXPLORATORY MODEL TESTING

##Load in, reference lists, and preprocessing

###Initialising reference objects for selecting relevant columns/indexes in functions
drops = ["AMP","SAM","TZP","CZO","CRO","CAZ","FEP","MEM","CIP","GEN","SXT","NIT","VAN" ]

###Load in
urines5 = pd.read_csv("urines5.csv")

###Ensure age is string
urines5['standard_age'] = urines5['standard_age'].map(str)

##Iterate over antibiotics and test all model types

###Initialise empty AUC df
big_auc_df = pd.DataFrame()

###Iterate
for antimicrobial in list(range(0, len(drops), 1)):

    ###et target ast outcome for antibiotic of choice
    y = urines5[drops[antimicrobial]]
    y = y.replace({'S': 1, 'R': 0}).astype(int)

    ###Dummy variables for predictor feature classes
    X = pd.get_dummies(urines5.drop(drops, axis=1))

    ###Train-test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=123)

    ###List of model types
    model_types = {
        'LogisticRegression': LogisticRegression(random_state=42, max_iter=1000),
        'SVC': SVC(kernel='rbf', probability=True, random_state=42),
        'KNeighborsClassifier': KNeighborsClassifier(n_neighbors=5),
        'DecisionTreeClassifier': DecisionTreeClassifier(random_state=42),
        'RandomForestClassifier': RandomForestClassifier(n_estimators=100, random_state=42),
        'GradientBoostingClassifier': GradientBoostingClassifier(random_state=42),
        'GaussianNB': GaussianNB(),
        'MLPClassifier': MLPClassifier(hidden_layer_sizes=(100,), max_iter=1000, random_state=42)
    }

    ###Initialise auc df
    auc_results = []

    ###Iterate over model types
    for name, clf in model_types.items():

        ###Fit model
        clf.fit(X_train, y_train)

        ###Output probability predictions
        y_pred_proba = clf.predict_proba(X_test)[:, 1]

        ###Get auc
        auroc = roc_auc_score(y_test, y_pred_proba)

        ###Add auc result to list
        auc_results.append({'Classifier': name, 'AUROC': auroc})

    ###Make dataframe of aucs for this antibiotic
    auc_results_df = pd.DataFrame(auc_results)
    auc_results_df['Antibiotic'] = drops[antimicrobial]

    ###Bind results to big auc df
    big_auc_df = pd.concat([big_auc_df, auc_results_df], ignore_index=True)

big_auc_df.to_csv('auc_testing_results.csv', index=False)
