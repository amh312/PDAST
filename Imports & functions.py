#IMPORTS AND FUNCTIONS

##Import from packages

###Data engineering
import pandas as pd
import numpy as np
from itertools import cycle
import pickle

###Scikit-learn
from sklearn.preprocessing import LabelBinarizer
from sklearn.metrics import auc,classification_report, roc_auc_score,roc_curve,RocCurveDisplay
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import cross_val_score, GridSearchCV, KFold

###Stats and data visualisation
import matplotlib.pyplot as plt
import statistics as stat
from mlxtend.evaluate import confusion_matrix

##Preprocessing

###Load-in and preprocessing
def df_clean(filepath,target,to_drop):
    # DATA LOAD-IN

    #read in urines df from csv
    df = pd.read_csv(filepath)

    # set target resistance var
    y = df[target]

    #make dummy vars for predictor vars
    df = pd.get_dummies(df.drop(to_drop, axis=1))

    #insert target var back into df
    df.insert(0, target, y)
    return(df)

##Hyperparameter tuning and feature selection

###LR hyperparameter tuning & feature selection
def lr_hyp_tune_feats(target,ht_features,test=0.2,random=1,folds=6,model_type="auto",cw=None,
                      scorer='f1_weighted',target_ab='antibiotic',targ_result='',analysis_type=''):

    #train-test split
    X_train, X_test, y_train, y_test = train_test_split(
        ht_features,
        target,
        test_size=test,
        random_state=random
    )

    #define log reg model
    log_reg = LogisticRegression(
        max_iter=2000,
        multi_class=model_type,
        class_weight=cw,
        penalty='l1',
        solver='liblinear'
    )

    #set 10-val c param grid
    param_grid = {"C": np.linspace(0.00001, 1, 10)}

    #set crossval parameters
    kf = KFold(n_splits=folds, random_state=random, shuffle=True)

    #define cv parameters
    logreg_cv = GridSearchCV(log_reg, param_grid, cv=kf,scoring=scorer)

    #run cv
    logreg_cv.fit(X_train, y_train)

    #print best c parameter value and score
    print("Tuned logreg parameters: {}".format(logreg_cv.best_params_))
    print("Tuned logreg score: {}".format(logreg_cv.best_score_))

    #report mean and sd of scores
    n_scores = cross_val_score(log_reg, X_train, y_train, scoring=scorer, cv=kf, n_jobs=-1)
    print('Mean Accuracy: %.3f (%.3f)' % (stat.mean(n_scores), stat.stdev(n_scores)))

    #retrieve best resistance model
    best_lr = logreg_cv.best_estimator_

    #coefs to glob envir
    global coefs

    #make df from coefs of best model and transpose
    coefs = pd.DataFrame(best_lr.coef_)
    coefs = coefs.transpose()

    #set ast result types as colnames if >1 result type
    if len(sorted(coefs)) > 1:
        coefs.columns = best_lr.classes_
    print("Total number of features:", len(ht_features.columns))

    #set var list to glob envir
    global vari_list
    vari_list = {}

    #check >1 ast result type for that antibiotic
    if len(sorted(coefs)) > 1:

        #populate predictor var list with nonzero predictors for each ast result type
        for label in best_lr.classes_:
            vari_list[label] = best_lr.feature_names_in_[coefs[coefs[label] != 0].index]
            print("Nonzero features for " + label + ": ", vari_list[label])
            print("Number of selected features for " + label + ": ",
                  np.count_nonzero(best_lr.feature_names_in_[coefs[coefs[label] != 0].index]))

            #set colours for plots depending on ast result
            if label == "R":
                line_col = "red"
            elif label == "I":
                line_col = "darkorange"
            elif label == "S":
                line_col = "green"
            elif label == "NT":
                line_col = "blue"

            #bar plot for positive coefficients for ast result
            if len(coefs.sort_values(by=0,ascending=False)[0:4] !=0):

                #set up for multiple bar plot for pred coeefs
                fig, ax = plt.subplots()

                #sort coefficients by value
                bar_coef = coefs.sort_values(by=0,ascending=False)[0:4]

                #for each ast result make a bar plot of predictor coefs
                bar_coef[label].plot.bar(color=line_col)

                #set feature names on axis and chart title
                ax.set_xticklabels(best_lr.feature_names_in_[coefs.sort_values(by=0,ascending=False)[0:4].index])
                ax.set_title(target_ab + ":\nStrongest positive coefficients for "+ label + " prediction")

                #set tight layout to ensure fit for pdf
                plt.tight_layout()

                #export to pdf and show figure
                plt.savefig(target_ab + label + "-posfeatures.pdf", format="pdf", bbox_inches="tight")
                plt.show()

            #same as above but for negative coefficients
            if len(coefs.sort_values(by=0,ascending=True)[0:4] != 0):

                #set for multiple plots
                fig, ax = plt.subplots()

                #set coefs all to positive and sort
                bar_coef = abs(coefs.sort_values(by=0,ascending=True)[0:4])

                #set colours according to ast result predicted
                bar_coef[label].plot.bar(color=line_col)

                #set axis labels and chart title
                ax.set_xticklabels(best_lr.feature_names_in_[coefs.sort_values(by=0,ascending=True)[0:4].index])
                ax.set_title(target_ab + ":\nStrongest negative coefficients for " + label + " prediction")

                #set tight layout and to pdf again
                plt.tight_layout()
                plt.savefig(target_ab + label + "-negfeatures.pdf", format="pdf", bbox_inches="tight")
                plt.show()

    #if there is only 1 ast result for that antibiotic
    else:

        #set vari list in same way as above
        vari_list = best_lr.feature_names_in_[(coefs[coefs != 0].dropna()).index]
        print("Nonzero features: ", vari_list)
        print("Number of selected features: ",
              np.count_nonzero(best_lr.feature_names_in_[(coefs[coefs != 0].dropna()).index]))

        #set space to have both plots
        fig, ax = plt.subplots()

        #set coefs for bar again, but this time drop missing vals
        bar_coef = coefs.sort_values(by=0,ascending=False)[0:4].dropna()

        #set plot as green for pos vals
        bar_coef.plot.bar(color="green")

        #set number of axis ticks manually and add labels, plus chart title
        plt.xticks(range(0, len(best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=False)[0:4].dropna()).index])),
                   best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=False)[0:4].dropna()).index])
        plt.title(target_ab + ":\nStrongest positive coefficients for susceptibility prediction")

        #legend
        plt.legend([])

        #tight layout again and pos features to bar plot pdf
        plt.tight_layout()
        plt.savefig(target_ab + targ_result + "-posfeatures.pdf", format="pdf", bbox_inches="tight")
        plt.show()

        #repeat for negative ceofficients
        fig, ax = plt.subplots()

        #set absolute vals of neg coefficients for bar purpose
        bar_coef = abs(coefs.sort_values(by=0,ascending=True)[0:4].dropna())

        #set neg features to red
        bar_coef.plot.bar(color="red")

        #manuallys set no of ticks and labels from predictor var names, plus title
        plt.xticks(range(0, len(best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=True)[0:4].dropna()).index])),
                   best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=True)[0:4].dropna()).index])
        plt.title(target_ab + ":\nStrongest negative coefficients for susceptibility prediction")

        #legend
        plt.legend([])

        #tight layout and to pdf
        plt.tight_layout()
        plt.savefig(target_ab + targ_result + "-negfeatures.pdf",format="pdf",bbox_inches="tight")
        plt.show()

    #set best c value in global envir
    global c_value

    #dictionary with best vals for params
    c_dict = logreg_cv.best_params_

    #get c value from dictionary
    c_value = list(c_dict.values())[0]

###LR hyperparameter tuning & feature selection (multinomial analysis)
def lr_hyp_tune_feats2(target,ht_features,test=0.2,random=1,folds=6,model_type="auto",cw=None,
                      scorer='f1_weighted',target_ab='antibiotic',targ_result=''):

    #train-test split
    X_train, X_test, y_train, y_test = train_test_split(
        ht_features,
        target,
        test_size=test,
        random_state=random
    )

    #set model parameters
    log_reg = LogisticRegression(
        max_iter=2000,
        multi_class=model_type,
        class_weight=cw,
        penalty='l1',
        solver='liblinear'
    )

    #set up values to test for c
    param_grid = {"C": np.linspace(0.00001, 1, 10)}

    #set up kfold params
    kf = KFold(n_splits=folds, random_state=random, shuffle=True)

    #set up grid search crossval
    logreg_cv = GridSearchCV(log_reg, param_grid, cv=kf,scoring=scorer)

    #run crossvalidation
    logreg_cv.fit(X_train, y_train)

    #report results of cv
    print("Tuned logreg parameters: {}".format(logreg_cv.best_params_))
    print("Tuned logreg score: {}".format(logreg_cv.best_score_))

    #get mean and sd of scores across crossval
    n_scores = cross_val_score(log_reg, X_train, y_train, scoring=scorer, cv=kf, n_jobs=-1)
    print('Mean Accuracy: %.3f (%.3f)' % (stat.mean(n_scores), stat.stdev(n_scores)))

    #get best model
    best_lr = logreg_cv.best_estimator_

    #set coefs to global envir
    global coefs

    #get coefs from best model as df and transpose
    coefs = pd.DataFrame(best_lr.coef_)
    coefs = coefs.transpose()

    #if >1 ast result for that ab then set ast results as coef df colnames
    if len(sorted(coefs)) > 1:
        coefs.columns = best_lr.classes_

    #print no. of features
    print("Total number of features:", len(ht_features.columns))

    #set variable list to global envir as dictionary
    global vari_list
    vari_list = {}

    #check >1 ast result again
    if len(sorted(coefs)) > 1:

        #iterate over ast results again
        for label in best_lr.classes_:

            #get coefficient variable names for each ast result type
            vari_list[label] = best_lr.feature_names_in_[coefs[coefs[label] != 0].index]
            print("Nonzero features for " + label + ": ", vari_list[label])
            print("Number of selected features for " + label + ": ",
                  np.count_nonzero(best_lr.feature_names_in_[coefs[coefs[label] != 0].index]))

            #set colours according to ast results
            if label == "R":
                line_col = "red"
            elif label == "I":
                line_col = "darkorange"
            elif label == "S":
                line_col = "green"
            elif label == "NT":
                line_col = "blue"

    #if only 1 ast result type
    else:

        #get coef names for non-zero predictor vars
        vari_list = best_lr.feature_names_in_[(coefs[coefs != 0].dropna()).index]
        print("Nonzero features: ", vari_list)
        print("Number of selected features: ",
              np.count_nonzero(best_lr.feature_names_in_[(coefs[coefs != 0].dropna()).index]))

        #set for multiple plots
        fig, ax = plt.subplots()

        #get top 5 non-zero coefficients
        bar_coef = coefs.sort_values(by=0,ascending=False)[0:4].dropna()

        #set green for pos coefficients
        bar_coef.plot.bar(color="green")

        #set x axis tick numbers and feature names for labels, and title
        plt.xticks(range(0, len(best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=False)[0:4].dropna()).index])),
                   best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=False)[0:4].dropna()).index])
        plt.title(target_ab + ":\nStrongest positive coefficients for susceptibility prediction")

        #legend
        plt.legend([])

        #tight layout and export to pdf / print plot
        plt.tight_layout()
        plt.savefig(target_ab + targ_result + "-posfeatures_multi.pdf", format="pdf", bbox_inches="tight")
        plt.show()

        #repeat for negative coefficients
        fig, ax = plt.subplots()

        #first 5 negatives
        bar_coef = abs(coefs.sort_values(by=0,ascending=True)[0:4].dropna())

        #set bars to red
        bar_coef.plot.bar(color="red")

        #x axis ticks and labels, title
        plt.xticks(range(0, len(best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=True)[0:4].dropna()).index])),
                   best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=True)[0:4].dropna()).index])
        plt.title(target_ab + ":\nStrongest negative coefficients for susceptibility prediction")

        #legens
        plt.legend([])

        #tight layout and export to pdf/print
        plt.tight_layout()
        plt.savefig(target_ab + targ_result + "-negfeatures_multi.pdf",format="pdf",bbox_inches="tight")
        plt.show()

    #c value in glob envir
    global c_value

    #get parameters from best model
    c_dict = logreg_cv.best_params_

    #get c value from start of best parameter list
    c_value = list(c_dict.values())[0]

###LR hyperparameter tuning & feature selection (out-of-sample time analysis)
def lr_hyp_tune_feats_t(target,ht_features,test=0.2,random=1,folds=6,model_type="auto",cw=None,
                      scorer='f1_weighted',target_ab='antibiotic',targ_result=''):

    #train_test split
    X_train, X_test, y_train, y_test = train_test_split(
        ht_features,
        target,
        test_size=test,
        random_state=random
    )

    #set log_reg model params
    log_reg = LogisticRegression(
        max_iter=2000,
        multi_class=model_type,
        class_weight=cw,
        penalty='l1',
        solver='liblinear'
    )

    #c values for tuning
    param_grid = {"C": np.linspace(0.00001, 1, 10)}

    #k fold params
    kf = KFold(n_splits=folds, random_state=random, shuffle=True)

    #gid cv params
    logreg_cv = GridSearchCV(log_reg, param_grid, cv=kf,scoring=scorer)

    #run cv and print results
    logreg_cv.fit(X_train, y_train)
    print("Tuned logreg parameters: {}".format(logreg_cv.best_params_))
    print("Tuned logreg score: {}".format(logreg_cv.best_score_))

    #get mean and sd of scores across cv
    n_scores = cross_val_score(log_reg, X_train, y_train, scoring=scorer, cv=kf, n_jobs=-1)
    print('Mean Accuracy: %.3f (%.3f)' % (stat.mean(n_scores), stat.stdev(n_scores)))

    #get best model
    best_lr = logreg_cv.best_estimator_

    #coefs to global as transposed df
    global coefs
    coefs = pd.DataFrame(best_lr.coef_)
    coefs = coefs.transpose()

    #set classes to colnames if enough ast results for df (i.e., at least 2)
    if len(sorted(coefs)) > 1:
        coefs.columns = best_lr.classes_
    print("Total number of features:", len(ht_features.columns))

    #global var list dictionary
    global vari_list
    vari_list = {}

    #check if >1 ast result
    if len(sorted(coefs)) > 1:

        #get non-zero feature names for each ast result prediction in vari liast
        for label in best_lr.classes_:
            vari_list[label] = best_lr.feature_names_in_[coefs[coefs[label] != 0].index]
            print("Nonzero features for " + label + ": ", vari_list[label])
            print("Number of selected features for " + label + ": ",
                  np.count_nonzero(best_lr.feature_names_in_[coefs[coefs[label] != 0].index]))

            #set plot colours according to AST results
            if label == "R":
                line_col = "red"
            elif label == "I":
                line_col = "darkorange"
            elif label == "S":
                line_col = "green"
            elif label == "NT":
                line_col = "blue"

            #check there are some non-zero coefficients
            if len(coefs.sort_values(by=0,ascending=False)[0:4] !=0):

                #set space for both bar plots
                fig, ax = plt.subplots()

                #set top 5 pos coefs as bar values
                bar_coef = coefs.sort_values(by=0,ascending=False)[0:4]

                #set colours of bars
                bar_coef[label].plot.bar(color="red")

                #set axis tick labels and title
                ax.set_xticklabels(best_lr.feature_names_in_[coefs.sort_values(by=0,ascending=False)[0:4].index])
                ax.set_title(target_ab + ":\nStrongest positive coefficients for "+ label + " prediction")

                #tight layout and out to pdf
                plt.tight_layout()
                plt.savefig(target_ab + label + "-posfeatures_time.pdf", format="pdf", bbox_inches="tight")


            if len(coefs.sort_values(by=0,ascending=True)[0:4] != 0):

                #multiple plot space
                fig, ax = plt.subplots()

                #absolutes of neg coefficients (top 5)
                bar_coef = abs(coefs.sort_values(by=0,ascending=True)[0:4])

                #pass coefs to bar plot
                bar_coef[label].plot.bar(color="red")

                #set axis tick labels and title
                ax.set_xticklabels(best_lr.feature_names_in_[coefs.sort_values(by=0,ascending=True)[0:4].index])
                ax.set_title(target_ab + ":\nStrongest negative coefficients for " + label + " prediction")

                #tight layout to fit titles and out to pdf
                plt.tight_layout()
                plt.savefig(target_ab + label + "-negfeatures_time.pdf", format="pdf", bbox_inches="tight")

    #if only 1 ast result
    else:

        #get nonzero coefs
        vari_list = best_lr.feature_names_in_[(coefs[coefs != 0].dropna()).index]
        print("Nonzero features: ", vari_list)
        print("Number of selected features: ",
              np.count_nonzero(best_lr.feature_names_in_[(coefs[coefs != 0].dropna()).index]))

        #multiple plot area
        fig, ax = plt.subplots()

        #get coef values
        bar_coef = coefs.sort_values(by=0,ascending=False)[0:4].dropna()

        #green bar colours for pos coefs
        bar_coef.plot.bar(color="green")

        #axis ticks and labs
        plt.xticks(range(0, len(best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=False)[0:4].dropna()).index])),
                   best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=False)[0:4].dropna()).index])
        plt.title(target_ab + ":\nStrongest positive coefficients for susceptibility prediction")

        #legend, tight layout and to pdf
        plt.legend([])
        plt.tight_layout()
        plt.savefig(target_ab + targ_result + "-posfeatures_time.pdf", format="pdf", bbox_inches="tight")

        #same for neg coefficients
        fig, ax = plt.subplots()

        #absolutes of negs
        bar_coef = abs(coefs.sort_values(by=0,ascending=True)[0:4].dropna())

        #red bars
        bar_coef.plot.bar(color="red")

        #axis ticks and labels
        plt.xticks(range(0, len(best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=True)[0:4].dropna()).index])),
                   best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=True)[0:4].dropna()).index])
        plt.title(target_ab + ":\nStrongest negative coefficients for susceptibility prediction")

        #set legend
        plt.legend([])

        #tight and to pdf, no print
        plt.tight_layout()
        plt.savefig(target_ab + targ_result + "-negfeatures_time.pdf",format="pdf",bbox_inches="tight")

    #c value to global
    global c_value

    #params from best model
    c_dict = logreg_cv.best_params_

    #pull out c hyperparam
    c_value = list(c_dict.values())[0]

##Model validation

###LR final validation run
def LR_multi_final(target,final_features,test=0.2,random=1,C_val=0.1,
             cw=None,model_type="auto",thr=0.5,
                   ab_of_interest="Ampicillin",av="micro",
                   analysis_type=''):

    #set things that will be needed outside function
    global log_reg
    global fairness_model
    global probs
    probs = {}
    global class_report
    global model_dict
    model_dict = {}
    global roc_auc
    global X_test

    #for multinomial analysis (i.e., s,r,i,nt)
    if len(target.unique()) > 2:

        #iterate over ast_results
        for label in vari_list:

            #train-test split
            X_train, X_test, y_train, y_test = train_test_split(
                final_features[vari_list[label]],
                target,
                test_size=test,
                random_state=random
            )

            #model training parameters
            log_reg = LogisticRegression(
                max_iter=2000,
                multi_class=model_type,
                C=C_val,
                class_weight=cw,
                penalty='l1',
                solver='liblinear'
            )

            #fit model
            log_reg.fit(X_train, y_train)

            #set model for later fairness analysis
            fairness_model = log_reg

            #put fit model in dictionary for later
            model_dict[label] = log_reg

            #get ast result classification predictions on test data
            y_pred = (log_reg.predict_proba(X_test)[:, 1] >= thr).astype(int)

            #same but for probability predictions
            y_pred_probs = log_reg.predict_proba(X_test)

            #get how many classes in data
            n_classes = len(np.unique(y_train))

            #put prob predictions into probs dictionary
            probs[label] = y_pred_probs

            #binarise labels to enable one-vs-rest approach
            label_binarizer = LabelBinarizer().fit(y_train)

            #binarised onehot encoding for test outcomes of interest
            y_onehot_test = label_binarizer.transform(y_test)

            #get dimensions to give no. of shapes and classes
            y_onehot_test.shape

            print('AUROC for '+label+': ',roc_auc_score(y_test, y_pred_probs,multi_class="ovr",average=av))

            #get class name of interest
            class_id = np.flatnonzero(label_binarizer.classes_ == label)[0]

            #set line colours according to ast results
            if label=="R":
                line_col="red"
            elif label=="I":
                line_col="darkorange"
            elif label=="S":
                line_col="green"
            elif label=="NT":
                line_col="blue"

            #make display for just result type of interest
            display = RocCurveDisplay.from_predictions(
                y_onehot_test[:, class_id],
                y_pred_probs[:, class_id],
                name=f"{label} vs the rest",
                color=line_col,
                plot_chance_level=True,
            )

            #set axis and roc plot titles
            _ = display.ax_.set(
                xlabel="False Positive Rate",
                ylabel="True Positive Rate",
                title="One-vs-Rest ROC curve for "+ab_of_interest+":\n"+label+" vs rest",
            )

            #save roc plot
            plt.savefig(ab_of_interest + label + analysis_type + "-roc.pdf", format="pdf", bbox_inches="tight")
            plt.show()

        #repeat train-test split
        X_train, X_test, y_train, y_test = train_test_split(
            final_features[vari_list["R"]],
            target,
            test_size=test,
            random_state=random
        )

        #reset model params
        log_reg = LogisticRegression(
            max_iter=2000,
            multi_class=model_type,
            C=C_val,
            class_weight=cw,
            penalty='l1',
            solver='liblinear'
        )

        #retrain model
        log_reg.fit(X_train, y_train)

        #set model to fairness
        fairness_model = log_reg

        #save resistance prediction model into dictionary
        model_dict["R"] = log_reg

        #get class predictions
        y_pred = (log_reg.predict_proba(X_test)[:, 1] >= thr).astype(int)

        #get prob predictions
        y_pred_probs = log_reg.predict_proba(X_test)

        #get number of ast result classes
        n_classes = len(np.unique(y_train))

        #get result classes themselves
        classes = log_reg.classes_

        #set resistance probability predictions in probs dictionary
        probs["R"] = y_pred_probs

        #get y predictions on test df
        y_testpred = log_reg.predict(X_test)

        #get classification report with full performance metrics
        class_report = classification_report(y_test, y_testpred,output_dict=True)

        #binarize labels for one-vs-rest analysis
        label_binarizer = LabelBinarizer().fit(y_train)

        #onehot encoding of outcome on test dataset
        y_onehot_test = label_binarizer.transform(y_test)

        #get no. of classes and cases
        y_onehot_test.shape

        #make micro-averaged roc-curve display across multinomial ast results
        display = RocCurveDisplay.from_predictions(
            y_onehot_test.ravel(),
            y_pred_probs.ravel(),
            name="micro-average OvR",
            color="purple",
            plot_chance_level=True,
        )

        #set title
        _ = display.ax_.set(
            xlabel="False Positive Rate",
            ylabel="True Positive Rate",
            title=ab_of_interest + " Micro-averaged One-vs-Rest\nROC",
        )

        #save pdf and print plot
        plt.savefig(ab_of_interest + analysis_type + "-avroc.pdf", format="pdf", bbox_inches="tight")
        plt.show()

        #store the fpr, tpr, and auroc for all averaging strategies
        fpr, tpr, roc_auc = dict(), dict(), dict()

        #output get micro-averaged roc curve and auroc
        fpr["micro"], tpr["micro"], _ = roc_curve(y_onehot_test.ravel(), y_pred_probs.ravel())
        roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

        print(f"Micro-averaged One-vs-Rest ROC AUC score:\n{roc_auc['micro']:.2f}")

        micro_roc_auc_ovr = roc_auc_score(
            y_test,
            y_pred_probs,
            multi_class="ovr",
            average="micro",
        )

        print(f"Micro-averaged One-vs-Rest ROC AUC score (2):\n{micro_roc_auc_ovr:.2f}")

        #macro-average
        for i in range(n_classes):

            #get auroc for each class
            fpr[i], tpr[i], _ = roc_curve(y_onehot_test[:, i], y_pred_probs[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])

        #vector of 1000 fpr values to plot across
        fpr_grid = np.linspace(0.0, 1.0, 1000)

        #matching array of zeroes to get tpr for each fpr
        mean_tpr = np.zeros_like(fpr_grid)

        #add x and y coordinates for roc plot to blank tpr array
        for i in range(n_classes):
            mean_tpr += np.interp(fpr_grid, fpr[i], tpr[i])

        #get means of tprs across ast result classes
        mean_tpr /= n_classes

        #set macro-averages in fpr and tpr dicts
        fpr["macro"] = fpr_grid
        tpr["macro"] = mean_tpr
        roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

        print(f"Macro-averaged One-vs-Rest ROC AUC score:\n{roc_auc['macro']:.2f}")


        #plot all roc curves together
        fig, ax = plt.subplots(figsize=(6, 6))

        #set plot for micro-averages (i.e., class-weighted means)
        plt.plot(
            fpr["micro"],
            tpr["micro"],
            label=f"micro-average ROC curve (AUC = {roc_auc['micro']:.2f})",
            color="deeppink",
            linestyle=":",
            linewidth=4,
        )

        #same for macr-averages (i.e., straight means)
        plt.plot(
            fpr["macro"],
            tpr["macro"],
            label=f"macro-average ROC curve (AUC = {roc_auc['macro']:.2f})",
            color="navy",
            linestyle=":",
            linewidth=4,
        )

        #set colour palette for four lines for result classes
        colors = cycle(["aqua", "darkorange", "cornflowerblue","limegreen"])

        #iterate across classes and respective colours to display rocs curves
        for class_id, color in zip(range(n_classes), colors):
            RocCurveDisplay.from_predictions(
                y_onehot_test[:, class_id],
                y_pred_probs[:, class_id],
                name=f"ROC curve for {classes[class_id]}",
                color=color,
                ax=ax,
                plot_chance_level=(class_id == 2),
            )

        #set axes and overall title
        _ = ax.set(
            xlabel="False Positive Rate",
            ylabel="True Positive Rate",
            title=ab_of_interest+" ROCs for\nOne-vs-Rest multiclass",
        )

        #save pdf
        plt.savefig(ab_of_interest + analysis_type + "-allroc.pdf", format="pdf", bbox_inches="tight")
        plt.show()

    #for main binomial analysis (i.e., s vs r)
    else:

        #train-test split
        X_train, X_test, y_train, y_test = train_test_split(
            final_features[vari_list],
            target,
            test_size=test,
            random_state=random
        )

        #training parameters
        log_reg = LogisticRegression(
            max_iter=2000,
            C=C_val,
            class_weight=cw,
            penalty='l1',
            solver='liblinear'
        )

        #fit model
        log_reg.fit(X_train, y_train)

        #assign model to model dict object
        model_dict = log_reg

        #assign fairness model for later
        fairness_model = log_reg

        #get classification predictions
        y_pred = (log_reg.predict_proba(X_test)[:, 0] >= thr).astype(int)

        #get probability predictions
        y_pred_probs = log_reg.predict_proba(X_test)[:,0]

        #assign probabilities to probs
        probs = y_pred_probs

        #get number of ast result classes
        n_classes = len(np.unique(y_train))

        #flip probs to ensure roc is predicting susceptibility side of binary prediction
        roc_auc = 1-roc_auc_score(y_test, y_pred_probs)

        #get outcome predictions on test dataset
        y_testpred = log_reg.predict(X_test)

        #get classification report for performance metrics
        class_report = classification_report(y_test, y_testpred,output_dict=True)

        #set line colours according to ast result types
        if log_reg.classes_[0] == "R":
            line_col = "red"
        elif log_reg.classes_[0] == "I":
            line_col = "darkorange"
        elif log_reg.classes_[0] == "S":
            line_col = "green"
        elif log_reg.classes_[0] == "NT":
            line_col = "blue"

        #set auroc display again but for binary prediction
        display = RocCurveDisplay.from_predictions(
            y_test,
            y_pred_probs,
            name=f"Binary prediction",
            color="purple",
            plot_chance_level=True,
            pos_label=log_reg.classes_[0]
        )

        #set axes names and title for binary prediction
        _ = display.ax_.set(
            xlabel="False Positive Rate",
            ylabel="True Positive Rate",
            title="Binary ROC curve for " + ab_of_interest + " susceptibility"
        )

        #export roc to pdf
        plt.savefig(ab_of_interest + analysis_type + "-roc.pdf", format="pdf", bbox_inches="tight")
        plt.show()

###LR final validation run (out-of-sample time analysis)
def LR_multi_final_t(target,final_features,test_datf,target2,test=0.2,random=1,C_val=0.1,
             cw=None,model_type="auto",thr=0.5,
                   ab_of_interest="Ampicillin",av="micro"):

    #set things that will be needed to global env
    global log_reg
    global fairness_model
    global probs
    probs = {}
    global class_report
    global model_dict
    model_dict = {}
    global roc_auc
    global X_test

    #for multinomial analysis (not actually used in out-of-sample time analysis)
    if len(target.unique()) > 2:

        #iterate over ast result types
        for label in vari_list:

            #train-tst split
            X_train, X_test, y_train, y_test = train_test_split(
                final_features[vari_list[label]],
                target,
                test_size=test,
                random_state=random
            )

            #model parameters
            log_reg = LogisticRegression(
                max_iter=2000,
                multi_class=model_type,
                C=C_val,
                class_weight=cw,
                penalty='l1',
                solver='liblinear'
            )

            #fit model
            log_reg.fit(X_train, y_train)

            #assign fairness for later
            fairness_model = log_reg

            #model for ast result of interest
            model_dict[label] = log_reg

            #classification predictions
            y_pred = (log_reg.predict_proba(X_test)[:, 1] >= thr).astype(int)

            #probability predictions
            y_pred_probs = log_reg.predict_proba(X_test)

            #number of ast result classes
            n_classes = len(np.unique(y_train))

            #prob predictions into probs dict
            probs[label] = y_pred_probs

            #binarise to enable one vs rest
            label_binarizer = LabelBinarizer().fit(y_train)

            #oneot encoding of binarised labels
            y_onehot_test = label_binarizer.transform(y_test)

            #dimensions for nuebr of cases and classes
            y_onehot_test.shape

            #return aroc
            print('AUROC score for '+label+': ',roc_auc_score(y_test, y_pred_probs,multi_class="ovr",average=av))

            #get indices of non-zero classes
            class_id = np.flatnonzero(label_binarizer.classes_ == label)[0]

            #set colours according to ast results predicted
            if label=="R":
                line_col="red"
            elif label=="I":
                line_col="darkorange"
            elif label=="S":
                line_col="green"
            elif label=="NT":
                line_col="blue"

            #display roc curve
            display = RocCurveDisplay.from_predictions(
                y_onehot_test[:, class_id],
                y_pred_probs[:, class_id],
                name=f"{label} vs the rest",
                color=line_col,
                plot_chance_level=True,
            )

            #set axes titles and overall title
            _ = display.ax_.set(
                xlabel="False Positive Rate",
                ylabel="True Positive Rate",
                title="One-vs-Rest ROC curve for "+ab_of_interest+":\n"+label+" vs rest",
            )

            #save auroc to pdf
            plt.savefig(ab_of_interest + label + "-roc_time.pdf", format="pdf", bbox_inches="tight")

        #train-test split again
        X_train, X_test, y_train, y_test = train_test_split(
            final_features[vari_list["R"]],
            target,
            test_size=test,
            random_state=random
        )

        #re-train model
        log_reg = LogisticRegression(
            max_iter=2000,
            multi_class=model_type,
            C=C_val,
            class_weight=cw,
            penalty='l1',
            solver='liblinear'
        )

        #fit model
        log_reg.fit(X_train, y_train)

        #fairness for later
        fairness_model = log_reg

        #set model to resistance prediction
        model_dict["R"] = log_reg

        #class predictions
        y_pred = (log_reg.predict_proba(X_test)[:, 1] >= thr).astype(int)

        #probability predictions
        y_pred_probs = log_reg.predict_proba(X_test)

        #number of ast result classes
        n_classes = len(np.unique(y_train))

        #list classes
        classes = log_reg.classes_

        #set prob predictions to probs df
        probs["R"] = y_pred_probs

        #get test predictions on test data
        y_testpred = log_reg.predict(X_test)

        #classificaiton report for performance metrics
        class_report = classification_report(y_test, y_testpred,output_dict=True)

        #generate one vs rest onehot labels with binarise again
        label_binarizer = LabelBinarizer().fit(y_train)
        y_onehot_test = label_binarizer.transform(y_test)
        y_onehot_test.shape  # (n_samples, n_classes)

        #make roc for micro-averaged results using ravel method
        display = RocCurveDisplay.from_predictions(
            y_onehot_test.ravel(),
            y_pred_probs.ravel(),
            name="micro-average OvR",
            color="purple",
            plot_chance_level=True,
        )

        #set axes labels and chart title
        _ = display.ax_.set(
            xlabel="False Positive Rate",
            ylabel="True Positive Rate",
            title=ab_of_interest + " Micro-averaged One-vs-Rest\nROC",
        )

        #save to pdf
        plt.savefig(ab_of_interest + "-avroc_time.pdf", format="pdf", bbox_inches="tight")


        #store fpr, tpr, and roc_auc for all averaging strategies as set of dicts
        fpr, tpr, roc_auc = dict(), dict(), dict()


        #compute micro-averaged roc and print
        fpr["micro"], tpr["micro"], _ = roc_curve(y_onehot_test.ravel(), y_pred_probs.ravel())
        roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

        print(f"Micro-averaged One-vs-Rest ROC AUC score:\n{roc_auc['micro']:.2f}")

        micro_roc_auc_ovr = roc_auc_score(
            y_test,
            y_pred_probs,
            multi_class="ovr",
            average="micro",
        )

        print(f"Micro-averaged One-vs-Rest ROC AUC score (2):\n{micro_roc_auc_ovr:.2f}")

        #repeat for macro-average
        for i in range(n_classes):
            fpr[i], tpr[i], _ = roc_curve(y_onehot_test[:, i], y_pred_probs[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])

        #fpr vals to plot across
        fpr_grid = np.linspace(0.0, 1.0, 1000)

        #set array of zeroes with fpr vector dimensions
        mean_tpr = np.zeros_like(fpr_grid)

        for i in range(n_classes):
            #set x and y coordinates for roc using interpolate
            mean_tpr += np.interp(fpr_grid, fpr[i], tpr[i])

        #get means (i.e., macro-average)
        mean_tpr /= n_classes

        #put values into dicts
        fpr["macro"] = fpr_grid
        tpr["macro"] = mean_tpr
        roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

        print(f"Macro-averaged One-vs-Rest ROC AUC score:\n{roc_auc['macro']:.2f}")


        #plot all roc curves together
        fig, ax = plt.subplots(figsize=(6, 6))

        #micro-averaged rocs
        plt.plot(
            fpr["micro"],
            tpr["micro"],
            label=f"micro-average ROC curve (AUC = {roc_auc['micro']:.2f})",
            color="deeppink",
            linestyle=":",
            linewidth=4,
        )

        #macro-averaged rocs
        plt.plot(
            fpr["macro"],
            tpr["macro"],
            label=f"macro-average ROC curve (AUC = {roc_auc['macro']:.2f})",
            color="navy",
            linestyle=":",
            linewidth=4,
        )

        #set colours
        colors = cycle(["aqua", "darkorange", "cornflowerblue","limegreen"])

        #iterate over classes and colours
        for class_id, color in zip(range(n_classes), colors):
            RocCurveDisplay.from_predictions(
                y_onehot_test[:, class_id],
                y_pred_probs[:, class_id],
                name=f"ROC curve for {classes[class_id]}",
                color=color,
                ax=ax,
                plot_chance_level=(class_id == 2),
            )

        #title and axes
        _ = ax.set(
            xlabel="False Positive Rate",
            ylabel="True Positive Rate",
            title=ab_of_interest+" ROCs for\nOne-vs-Rest multiclass",
        )

        plt.savefig(ab_of_interest + "-allroc_time.pdf", format="pdf", bbox_inches="tight")

    #binomial lr (i.e., s vs r) - used in actual out-of-sample time analysis
    else:

        #train-test split to get training dataset from time period 1
        X_train, no_X_test, y_train, no_y_test = train_test_split(
            final_features[vari_list],
            target,
            test_size=test,
            random_state=random
        )

        #train-test split  to get testing dataset for time period 2
        X_no_train, X_test, y_no_train, y_test = train_test_split(
            test_datf[vari_list],
            target2,
            test_size=test,
            random_state=random
        )

        #set model parameters
        log_reg = LogisticRegression(
            max_iter=2000,
            C=C_val,
            class_weight=cw,
            penalty='l1',
            solver='liblinear'
        )

        #fit model
        log_reg.fit(X_train, y_train)

        #set to model_dict and set fairness model for later
        model_dict = log_reg
        fairness_model = log_reg

        #get class predictions
        y_pred = (log_reg.predict_proba(X_test)[:, 0] >= thr).astype(int)

        #get probability predictions and set to probs object
        y_pred_probs = log_reg.predict_proba(X_test)[:,0]
        probs = y_pred_probs

        #get n classes
        n_classes = len(np.unique(y_train))
        print('ROC_AUC score: ', 1-roc_auc_score(y_test, y_pred_probs))

        #flip roc to ensure predicting susceptibility
        roc_auc = 1-roc_auc_score(y_test, y_pred_probs)

        #class predictions for class report
        y_testpred = log_reg.predict(X_test)

        #classification report for performance metrics
        class_report = classification_report(y_test, y_testpred,output_dict=True)

        #set line colours as per ast result type
        if log_reg.classes_[0] == "R":
            line_col = "red"
        elif log_reg.classes_[0] == "I":
            line_col = "darkorange"
        elif log_reg.classes_[0] == "S":
            line_col = "green"
        elif log_reg.classes_[0] == "NT":
            line_col = "blue"

        #roc curve for binary prediction on time period train-test pair of interest
        display = RocCurveDisplay.from_predictions(
            y_test,
            y_pred_probs,
            name=f"Binary prediction",
            color="purple",
            plot_chance_level=True,
            pos_label=log_reg.classes_[0]
        )

        #axes and title
        _ = display.ax_.set(
            xlabel="False Positive Rate",
            ylabel="True Positive Rate",
            title="Binary ROC curve for " + ab_of_interest + " susceptibility"
        )

        #to pdf
        plt.savefig(ab_of_interest + "-roc_time.pdf", format="pdf", bbox_inches="tight")

###Compiling performance metrics
def result_compiler(data,output_filename):

    #make empty list
    flat_data = []

    #iterate over model and metrics to add performance metrics for each model
    for model, metrics in data.items():
        for metric, values in metrics.items():

            #if values are a dictionary then do additional iteration loop across sub-metrics
            if isinstance(values, dict):
                for sub_metric, score in values.items():
                    #append models and metrics to list
                    flat_data.append([model, metric, sub_metric, score])

            #if values are a vector then append straight to list
            else:
                flat_data.append([model, metric, None, values])

    #convert list to df
    df = pd.DataFrame(flat_data, columns=['Model', 'Metric', 'Sub-metric', 'Score'])

    #susceptible result performance metrics
    s_df = df[(df['Metric'] == 'S')]
    s_df = s_df.reset_index(drop=True)

    #resistant result performance metrics
    r_df = df[(df['Metric'] == 'R')]
    r_df = r_df.reset_index(drop=True)

    #accuracy performance metrics
    acc_df = df[(df['Metric'] == 'accuracy')]
    acc_df = acc_df.reset_index(drop=True)

    #bind s,r,acc vectors into single df
    df = pd.concat([s_df,r_df,acc_df])

    #export performance metrics df to csv
    df.to_csv(output_filename, index=False)

    return(df)

##Model fairness assessment

###Performing analysis
def mod_fairness(protvarv,chosen_model,target,final_features,test=0.2,random=1,
                 thr=0.5,av="micro",Bool=True):

    global log_reg
    global probs
    probs = {}
    global class_report
    global model_dict
    model_dict = {}

    global X_test
    global y_test

    if len(target.unique()) > 2:

        for label in vari_list:

            X_train, X_test, y_train, y_test = train_test_split(
                final_features,
                target,
                test_size=test,
                random_state=random
            )

            X_test = X_test[X_test[protvarv] == Bool]
            y_test = y_test[X_test[X_test[protvarv] == Bool].index]
            X_test = X_test[vari_list[label]]


            # Train a new model and evaluate performance

            log_reg = chosen_model
            y_pred = (log_reg.predict_proba(X_test)[:, 1] >= thr).astype(int)
            y_pred_probs = log_reg.predict_proba(X_test)
            n_classes = len(np.unique(y_train))
            probs[label] = y_pred_probs

            # Generate ROC curve values: fpr, tpr, thresholds
            label_binarizer = LabelBinarizer().fit(y_train)
            y_onehot_test = label_binarizer.transform(y_test)
            y_onehot_test.shape  # (n_samples, n_classes)

            print('ROC_AUC score for '+label+': ',roc_auc_score(y_test, y_pred_probs,multi_class="ovr",average=av))

            class_id = np.flatnonzero(label_binarizer.classes_ == label)[0]

            if label=="R":
                line_col="red"
            elif label=="I":
                line_col="darkorange"
            elif label=="S":
                line_col="green"
            elif label=="NT":
                line_col="blue"


        X_train, X_test, y_train, y_test = train_test_split(
            final_features,
            target,
            test_size=test,
            random_state=random
        )

        X_test = X_test[X_test[protvarv] == Bool]
        y_test = y_test[X_test[X_test[protvarv] == Bool].index]
        X_test = X_test[vari_list['R']]

        # Train a new model and evaluate performance

        log_reg = chosen_model
        y_pred = (log_reg.predict_proba(X_test)[:, 1] >= thr).astype(int)
        y_pred_probs = log_reg.predict_proba(X_test)
        n_classes = len(np.unique(y_train))
        classes = log_reg.classes_
        probs["R"] = y_pred_probs
        y_testpred = log_reg.predict(X_test)
        cm = confusion_matrix(y_test, y_testpred)
        TN, FP, FN, TP = cm.ravel()
        N = TP + FP + FN + TN  # Total population
        ACC = (TP + TN) / N  # Accuracy
        TPR = TP / (TP + FN)  # True positive rate
        FPR = FP / (FP + TN)  # False positive rate
        FNR = FN / (TP + FN)  # False negative rate
        PPP = (TP + FP) / N  # % predicted as positive

        class_report = np.array([ACC, TPR, FPR, FNR, PPP])

        # Generate ROC curve values: fpr, tpr, thresholds
        label_binarizer = LabelBinarizer().fit(y_train)
        y_onehot_test = label_binarizer.transform(y_test)
        y_onehot_test.shape  # (n_samples, n_classes)

        # store the fpr, tpr, and roc_auc for all averaging strategies
        fpr, tpr, roc_auc = dict(), dict(), dict()


        # Compute micro-average ROC curve and ROC area
        fpr["micro"], tpr["micro"], _ = roc_curve(y_onehot_test.ravel(), y_pred_probs.ravel())
        roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

        print(f"Micro-averaged One-vs-Rest ROC AUC score:\n{roc_auc['micro']:.2f}")

        micro_roc_auc_ovr = roc_auc_score(
            y_test,
            y_pred_probs,
            multi_class="ovr",
            average="micro",
        )

        print(f"Micro-averaged One-vs-Rest ROC AUC score (2):\n{micro_roc_auc_ovr:.2f}")

        #macro-average

        for i in range(n_classes):
            fpr[i], tpr[i], _ = roc_curve(y_onehot_test[:, i], y_pred_probs[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])

        fpr_grid = np.linspace(0.0, 1.0, 1000)

        # Interpolate all ROC curves at these points
        mean_tpr = np.zeros_like(fpr_grid)

        for i in range(n_classes):
            mean_tpr += np.interp(fpr_grid, fpr[i], tpr[i])  # linear interpolation

        # Average it and compute AUC
        mean_tpr /= n_classes

        fpr["macro"] = fpr_grid
        tpr["macro"] = mean_tpr
        roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

        print(f"Macro-averaged One-vs-Rest ROC AUC score:\n{roc_auc['macro']:.2f}")

    else:


        X_train, X_test, y_train, y_test = train_test_split(
            final_features,
            target,
            test_size=test,
            random_state=random
        )

        X_test = X_test[X_test[protvarv] == Bool]
        y_test = y_test[X_test[X_test[protvarv] == Bool].index]
        X_test = X_test[vari_list]

        # Train a new model and evaluate performance

        log_reg = chosen_model

        y_pred = (log_reg.predict_proba(X_test)[:, 0] >= thr).astype(int)
        y_pred_probs = log_reg.predict_proba(X_test)[:,0]
        probs = y_pred_probs
        n_classes = len(np.unique(y_train))
        y_testpred = log_reg.predict(X_test)
        cm = confusion_matrix(y_test, y_testpred)
        TN, FP, FN, TP = cm.ravel()
        N = TP + FP + FN + TN  # Total population
        ACC = (TP + TN) / N  # Accuracy
        TPR = TP / (TP + FN)  # True positive rate
        FPR = FP / (FP + TN)  # False positive rate
        FNR = FN / (TP + FN)  # False negative rate
        PPP = (TP + FP) / N  # % predicted as positive

        class_report = np.array([ACC, TPR, FPR, FNR, PPP])

        if log_reg.classes_[0] == "R":
            line_col = "red"
        elif log_reg.classes_[0] == "I":
            line_col = "darkorange"
        elif log_reg.classes_[0] == "S":
            line_col = "green"
        elif log_reg.classes_[0] == "NT":
            line_col = "blue"

###Extracting fairness metrics
def fairness_printer(df, metric_num):
    results = {}
    for prot in prot_vars:
        if metric_num < len(df[prot]):
            results[prot] = round(df[prot][metric_num], 2)
        else:
            results[prot] = None

    return results

###Compiling fairness metrics
def compile_fairness_results(fairness_df,output_filename):

    metrictable = {}
    for antimicrobial in drops:
        metrictable[antimicrobial] = {}
        for metricnumber, metricname in enumerate(metrics_fairness):
            metrictable[antimicrobial][metricname] = fairness_printer(fairness_df[antimicrobial], metricnumber)

    df_list = []

    for antimicrobial in metrictable:
        for metric in metrictable[antimicrobial]:
            row = {
                'Antimicrobial': antimicrobial,
                'Metric': metric
            }
            row.update(metrictable[antimicrobial][metric])
            df_list.append(row)

    result_df = pd.DataFrame(df_list)

    result_df.to_csv(output_filename, index=False)

    return(result_df)

##Stability sensitivity analysis

###Model stability assessment
def mod_stab_pred(target,final_features,test=0.2,random=1,C_val=0.1,
             cw=None,model_type="auto",thr=0.5,
                   abx="Ampicillin",av="micro",to_drop=""):

    global size_probs
    size_probs = {}

    global size_aucs
    size_aucs ={}
    global classes
    classes = {}

    test_sizes = [0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98]

    for test_sizer in test_sizes:

        if len(target.unique()) > 2:

            X_train, X_test, y_train, y_test = train_test_split(
                final_features[vari_list['R']],
                target,
                test_size=test_sizer,
                random_state=random
            )

            breaker = 1

            while len(y_train.unique()) < 2 or len(y_test.unique()) < 2:
                print("only one class. trying again...")

                X_train, X_test, y_train, y_test = train_test_split(
                    final_features[vari_list['R']],
                    target,
                    test_size=test_sizer,
                    random_state=random + 1
                )
                breaker = breaker + 1
                if breaker == 10:
                    print("giving up.")
                    break

            # Train a new model and evaluate performance

            log_reg = LogisticRegression(
                max_iter=2000,
                multi_class=model_type,
                C=C_val,
                class_weight=cw,
                penalty = 'l1',
                solver = 'liblinear'
            )
            log_reg.fit(X_train, y_train)
            classes[str(test_sizer)] = log_reg.classes_
            y_pred_probs = log_reg.predict_proba(X_test)
            size_probs[str(test_sizer)] = y_pred_probs
            size_probs['Antimicrobial'] = abx
            test_extra = 'actval' + str(test_sizer)
            size_probs[test_extra] = np.array(y_test)
            if len(y_test.unique()) == list(y_pred_probs.shape)[1]:
                size_aucs[test_sizer] = roc_auc_score(y_test, y_pred_probs,multi_class="ovr",average=av)
                size_aucs['Antimicrobial'] = abx
            else:
                size_aucs[test_sizer] = float("NaN")
                size_aucs['Antimicrobial'] = abx

        else:

            X_train, X_test, y_train, y_test = train_test_split(
                final_features[vari_list],
                target,
                test_size=test_sizer,
                random_state=random
            )

            breaker = 1

            while len(y_train.unique()) < 2 or len(y_test.unique()) < 2:
                print("only one class. trying again...")

                X_train, X_test, y_train, y_test = train_test_split(
                    final_features[vari_list],
                    target,
                    test_size=test_sizer,
                    random_state=random+1
                )
                breaker = breaker+1
                if breaker ==10:
                    print("giving up.")
                    break


            # Train a new model and evaluate performance

            log_reg = LogisticRegression(
                max_iter=2000,
                C=C_val,
                class_weight=cw,
                penalty='l1',
                solver='liblinear'
            )
            log_reg.fit(X_train, y_train)
            y_pred_probs = log_reg.predict_proba(X_test)
            classes[str(test_sizer)] = log_reg.classes_
            size_probs[str(test_sizer)] = y_pred_probs
            size_probs['Antimicrobial'] = abx
            test_extra = 'actval' + str(test_sizer)
            size_probs[test_extra] = np.array(y_test)

            if len(y_test.unique()) == list(y_pred_probs.shape)[1]:
                size_aucs[test_sizer] = roc_auc_score(y_test, y_pred_probs[:,0])
                size_aucs['Antimicrobial'] = abx
            else:
                size_aucs[test_sizer] = float("NaN")
                size_aucs['Antimicrobial'] = abx

##Microsimulation study

###Probability predictions
def prob_predict(abx,abx_name,to_drop,test_filepath,var_filepath,fit_filepath,
                 full_varlist,result_optimised_for,ref):
    global test_df
    global vari_list
    global model_dict
    global X_test

    test_df = pd.read_csv(test_filepath)
    test_ids = test_df[to_drop]
    test_df = pd.get_dummies(test_df.drop(to_drop, axis=1))
    test_df = test_df.reindex(columns=full_varlist)
    test_df = test_df.fillna(False)

    with open(var_filepath, 'rb') as f:
        vari_list = pickle.load(f)
    with open(fit_filepath, 'rb') as f:
        model_dict = pickle.load(f)

    global probs
    probs = {}

    global probs_df

    global probs_df_final

    if len(ref[abx_name].unique()) > 2:

        for label in vari_list:
            X_test = test_df[vari_list[label]]
            y_pred_probs = model_dict[label].predict_proba(X_test)
            probs[label] = y_pred_probs

        probs_df = pd.DataFrame()
        for label in probs:
            df2 = pd.DataFrame.from_dict(probs[label])
            df2['Variable_tuning'] = label
            probs_df = pd.concat([probs_df, df2])

        names = sorted(probs)
        names.append("Variable_tuning")
        probs_df.columns = names
        probs_df.to_csv("probs_df.csv")

        probs_df_final = pd.DataFrame()

        df_filter = result_optimised_for
        probs_df_final = probs_df.query('Variable_tuning==@df_filter')

        probs_df_final = pd.concat([test_ids, probs_df_final], axis=1)
        probs_df_final['Antimicrobial'] = abx

    else:

        X_test = test_df[vari_list]
        y_pred_probs = model_dict.predict_proba(X_test)
        probs = y_pred_probs

        probs_df = pd.DataFrame.from_dict(probs)

        probs_df.to_csv("probs_df.csv")

        probs_df_final = probs_df
        probs_df_final.columns = model_dict.classes_

        probs_df_final["Variable tuning"] = model_dict.classes_[1]

        probs_df_final = pd.concat([test_ids, probs_df_final], axis=1)
        probs_df_final['Antimicrobial'] = abx

##Out-of-sample time analysis

###Within own time period
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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

###Across time periods
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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
        lr_hyp_tune_feats_t(target=urines5["Y"],
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

##Wrapping functions

###Primary test/train/validate
def train_test_validate_main(csv,ab,abx,antibiotic,type=''):

    global features
    global urines5

    urines5 = pd.read_csv(csv)
    urines5['standard_age'] = urines5['standard_age'].map(str)
    y = urines5[abx]
    urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
    urines5.insert(0, abx, y)


    #set Y variable
    urines5['Y'] = urines5[abx]
    urines5 = urines5.drop(abx,axis=1)
    features = urines5.drop('Y',axis=1)

    #C hyperparameter tuning and feature selection
    lr_hyp_tune_feats(target=urines5["Y"],
                ht_features=features,
                test=0.2,
                random=100,
                folds=6,
                model_type="ovr",
                cw='balanced',
                scorer='roc_auc_ovr',
                target_ab = antibiotic,
                      analysis_type=type)

    vari_overall[abx] = vari_list
    c_values[abx] = c_value

    #LR final run
    LR_multi_final(target=urines5["Y"],
             final_features=features,
             test=0.2,
             random=100,
             C_val=c_value,
             cw='balanced',
             model_type="ovr",
             thr=0.5,
             ab_of_interest=antibiotic,
             av="micro",
                   analysis_type=type)

    class_reps[abx] = class_report

    #save fit, variables, and hyperparameters
    with open(ab+type+'file.pickle','wb') as f:
        pickle.dump(vari_list,f)
    with open(ab+type+'fits.pickle','wb') as f:
        pickle.dump(model_dict,f)

###Multinomial test/train/validate
def train_test_validate_multinomial(csv,ab,abx,antibiotic):

    global features
    global urines5

    urines5 = pd.read_csv(csv)
    urines5['standard_age'] = urines5['standard_age'].map(str)
    y = urines5[abx]
    urines5 = pd.get_dummies(urines5.drop(drops, axis=1))
    urines5.insert(0, abx, y)


    #set Y variable
    urines5['Y'] = urines5[abx]
    urines5 = urines5.drop(abx,axis=1)
    features = urines5.drop('Y',axis=1)


    #C hyperparameter tuning and feature selection
    lr_hyp_tune_feats2(target=urines5["Y"],
                ht_features=features,
                test=0.2,
                random=100,
                folds=6,
                model_type="ovr",
                cw='balanced',
                scorer='roc_auc_ovr',
                target_ab = antibiotic)

    vari_overall[abx] = vari_list
    c_values[abx] = c_value

    #LR final run
    LR_multi_final(target=urines5["Y"],
             final_features=features,
             test=0.2,
             random=100,
             C_val=c_value,
             cw='balanced',
             model_type="ovr",
             thr=0.5,
             ab_of_interest=antibiotic,
             av="micro",
                   analysis_type="multinomial")

    class_reps[abx] = class_report

    #save fit, variables, and hyperparameters
    with open(ab+'file_multi.pickle','wb') as f:
        pickle.dump(vari_list,f)
    with open(ab+'fits_multi.pickle','wb') as f:
        pickle.dump(model_dict,f)

###Model fairness analysis
def iterate_mod_fairness(df,feature_df, antibiotic,abx,ab):

    global overallfairT
    global overallfairF

    ####Initialising reference objects for fairness analysis
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

    fair_classrepsT = {}
    fair_classrepsF = {}

    for protvar in prot_vars:
        mod_fairness(protvarv=protvar,
                     chosen_model=fairness_model,
                     target=df["Y"],
                     final_features=feature_df,
                     test=0.2,
                     random=100,
                     thr=0.5,
                     av="micro",
                     Bool=True)

        fair_classrepsT[protvar] = class_report

        mod_fairness(protvarv=protvar,
                     chosen_model=fairness_model,
                     target=df["Y"],
                     final_features=feature_df,
                     test=0.2,
                     random=100,
                     thr=0.5,
                     av="micro",
                     Bool=False)

        fair_classrepsF[protvar] = class_report

    with open(abx + 'fairF.pickle', 'wb') as f:
        pickle.dump(fair_classrepsF, f)
    with open(abx + 'fairT.pickle', 'wb') as f:
        pickle.dump(fair_classrepsT, f)

    overallfairT[ab] = fair_classrepsT
    overallfairF[ab] = fair_classrepsF

###Stability assessment
def iter_stab_pred(df,feature_df,Antimicrobial_agent,abx):

        global p1
        global p2
        global p3
        global p4
        global p5
        global p6
        global p7
        global p8
        global p9

        p2 = pd.DataFrame()
        p3 = pd.DataFrame()
        p4 = pd.DataFrame()
        p5 = pd.DataFrame()
        p6 = pd.DataFrame()
        p7 = pd.DataFrame()
        p8 = pd.DataFrame()
        p9 = pd.DataFrame()

        for i in list(range(1, 101, 1)):
            mod_stab_pred(target=df["Y"],
                          final_features=feature_df,
                          test=0.2,
                          random=i,
                          C_val=c_value,
                          cw='balanced',
                          model_type="ovr",
                          thr=0.5,
                          abx=Antimicrobial_agent,
                          av="micro",
                          to_drop=to_drop)

            pointtwo = pd.DataFrame.from_dict(size_probs['0.84'])
            pointtwo.columns = classes['0.84']
            pointtwob = pd.DataFrame.from_dict(size_probs['actval0.84'])
            pointtwo['actval'] = pointtwob
            pointtwo['AUC'] = size_aucs[0.84]
            pointtwo['Antimicrobial'] = Antimicrobial_agent
            pointtwo['model'] = i
            p2 = pd.concat([p2, pointtwo])

            pointthree = pd.DataFrame.from_dict(size_probs['0.86'])
            pointthree.columns = classes['0.86']
            pointthreeb = pd.DataFrame.from_dict(size_probs['actval0.86'])
            pointthree['actval'] = pointthreeb
            pointthree['AUC'] = size_aucs[0.86]
            pointthree['Antimicrobial'] = Antimicrobial_agent
            pointthree['model'] = i
            p3 = pd.concat([p3, pointthree])

            pointfour = pd.DataFrame.from_dict(size_probs['0.88'])
            pointfour.columns = classes['0.88']
            pointfourb = pd.DataFrame.from_dict(size_probs['actval0.88'])
            pointfour['actval'] = pointfourb
            pointfour['AUC'] = size_aucs[0.88]
            pointfour['Antimicrobial'] = Antimicrobial_agent
            pointfour['model'] = i
            p4 = pd.concat([p4, pointfour])

            pointfive = pd.DataFrame.from_dict(size_probs['0.9'])
            pointfive.columns = classes['0.9']
            pointfiveb = pd.DataFrame.from_dict(size_probs['actval0.9'])
            pointfive['actval'] = pointfiveb
            pointfive['AUC'] = size_aucs[0.9]
            pointfive['Antimicrobial'] = Antimicrobial_agent
            pointfive['model'] = i
            p5 = pd.concat([p5, pointfive])

            pointsix = pd.DataFrame.from_dict(size_probs['0.92'])
            pointsix.columns = classes['0.92']
            pointsixb = pd.DataFrame.from_dict(size_probs['actval0.92'])
            pointsix['actval'] = pointsixb
            pointsix['AUC'] = size_aucs[0.92]
            pointsix['Antimicrobial'] = Antimicrobial_agent
            pointsix['model'] = i
            p6 = pd.concat([p6, pointsix])

            pointseven = pd.DataFrame.from_dict(size_probs['0.94'])
            pointseven.columns = classes['0.94']
            pointsevenb = pd.DataFrame.from_dict(size_probs['actval0.94'])
            pointseven['actval'] = pointsevenb
            pointseven['AUC'] = size_aucs[0.94]
            pointseven['Antimicrobial'] = Antimicrobial_agent
            pointseven['model'] = i
            p7 = pd.concat([p7, pointseven])

            pointeight = pd.DataFrame.from_dict(size_probs['0.96'])
            pointeight.columns = classes['0.96']
            pointeightb = pd.DataFrame.from_dict(size_probs['actval0.96'])
            pointeight['actval'] = pointeightb
            pointeight['AUC'] = size_aucs[0.96]
            pointeight['Antimicrobial'] = Antimicrobial_agent
            pointeight['model'] = i
            p8 = pd.concat([p8, pointeight])

            pointnine = pd.DataFrame.from_dict(size_probs['0.98'])
            pointnine.columns = classes['0.98']
            pointnineb = pd.DataFrame.from_dict(size_probs['actval0.98'])
            pointnine['actval'] = pointnineb
            pointnine['AUC'] = size_aucs[0.98]
            pointnine['Antimicrobial'] = Antimicrobial_agent
            pointnine['model'] = i
            p9 = pd.concat([p9, pointnine])

            p2.to_csv("p2_"+abx+".csv")
            p3.to_csv("p3_"+abx+".csv")
            p4.to_csv("p4_"+abx+".csv")
            p5.to_csv("p5_"+abx+".csv")
            p6.to_csv("p6_"+abx+".csv")
            p7.to_csv("p7_"+abx+".csv")
            p8.to_csv("p8_"+abx+".csv")
            p9.to_csv("p9_"+abx+".csv")

            print(i)
