#IMPORT
import random
from itertools import cycle
import matplotlib.pyplot as plt
from matplotlib import transforms,pyplot
import pandas as pd
import numpy as np
from itertools import compress
from statistics import mean, median,mode
import sklearn.metrics
from sklearn.metrics import mean_absolute_error
from sklearn.preprocessing import LabelBinarizer
from sklearn.preprocessing import LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFE
from mlxtend.feature_selection import sequential_feature_selector
from sklearn.ensemble import VotingClassifier, BaggingClassifier, GradientBoostingClassifier, RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.model_selection import train_test_split, RepeatedStratifiedKFold, cross_val_score
import seaborn as sns
from sklearn.metrics import auc, f1_score, precision_score, recall_score, accuracy_score, balanced_accuracy_score,classification_report, roc_auc_score,roc_curve,make_scorer,RocCurveDisplay
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score, GridSearchCV, KFold, RandomizedSearchCV
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
import statistics as stat
from mlxtend.evaluate import confusion_matrix
from mlxtend.plotting import plot_confusion_matrix
from statistics import stdev,mean
import pickle
from joblib import dump,load

# probability prediction
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

#data load in, check and cleaning
def df_clean(filepath,target,to_drop):
    # DATA LOAD-IN

    df = pd.read_csv(filepath)
    # convert to dummy variables
    y = df[target]
    df = pd.get_dummies(df.drop(to_drop, axis=1))
    df.insert(0, target, y)
    return(df)



#LR hyperparameter tuning & feature selection
def lr_hyp_tune_feats(df,target,ht_features,test=0.2,random=1,folds=6,model_type="auto",cw=None,
                      scorer='f1_weighted',target_ab='antibiotic',targ_result=''):
    X_train, X_test, y_train, y_test = train_test_split(
        ht_features,
        target,
        test_size=test,
        random_state=random
    )
    log_reg = LogisticRegression(
        max_iter=2000,
        multi_class=model_type,
        class_weight=cw,
        penalty='l1',
        solver='liblinear'
    )
    param_grid = {"C": np.linspace(0.00001, 1, 10)}
    kf = KFold(n_splits=folds, random_state=random, shuffle=True)
    logreg_cv = GridSearchCV(log_reg, param_grid, cv=kf,scoring=scorer)
    logreg_cv.fit(X_train, y_train)
    print("Tuned logreg parameters: {}".format(logreg_cv.best_params_))
    print("Tuned logreg score: {}".format(logreg_cv.best_score_))
    n_scores = cross_val_score(log_reg, X_train, y_train, scoring=scorer, cv=kf, n_jobs=-1)
    print('Mean Accuracy: %.3f (%.3f)' % (stat.mean(n_scores), stat.stdev(n_scores)))
    best_lr = logreg_cv.best_estimator_
    global coefs
    coefs = pd.DataFrame(best_lr.coef_)
    coefs = coefs.transpose()
    if len(sorted(coefs)) > 1:
        coefs.columns = best_lr.classes_
    print("Total number of features:", len(ht_features.columns))
    global vari_list
    vari_list = {}

    if len(sorted(coefs)) > 1:
        for label in best_lr.classes_:
            vari_list[label] = best_lr.feature_names_in_[coefs[coefs[label] != 0].index]
            print("Nonzero features for " + label + ": ", vari_list[label])
            print("Number of selected features for " + label + ": ",
                  np.count_nonzero(best_lr.feature_names_in_[coefs[coefs[label] != 0].index]))

            if label == "R":
                line_col = "red"
            elif label == "I":
                line_col = "darkorange"
            elif label == "S":
                line_col = "green"
            elif label == "NT":
                line_col = "blue"

            if len(coefs.sort_values(by=0,ascending=False)[0:4] !=0):

                fig, ax = plt.subplots()
                bar_coef = coefs.sort_values(by=0,ascending=False)[0:4]
                bar_coef[label].plot.bar(color=line_col)
                ax.set_xticklabels(best_lr.feature_names_in_[coefs.sort_values(by=0,ascending=False)[0:4].index])
                ax.set_title(target_ab + ":\nStrongest positive coefficients for "+ label + " prediction")
                plt.tight_layout()
                plt.savefig(target_ab + label + "-posfeatures.pdf", format="pdf", bbox_inches="tight")
                plt.show()

            if len(coefs.sort_values(by=0,ascending=True)[0:4] != 0):

                fig, ax = plt.subplots()
                bar_coef = abs(coefs.sort_values(by=0,ascending=True)[0:4])
                bar_coef[label].plot.bar(color=line_col)
                ax.set_xticklabels(best_lr.feature_names_in_[coefs.sort_values(by=0,ascending=True)[0:4].index])
                ax.set_title(target_ab + ":\nStrongest negative coefficients for " + label + " prediction")
                plt.tight_layout()
                plt.savefig(target_ab + label + "-negfeatures.pdf", format="pdf", bbox_inches="tight")
                plt.show()

    else:
        vari_list = best_lr.feature_names_in_[(coefs[coefs != 0].dropna()).index]
        print("Nonzero features: ", vari_list)
        print("Number of selected features: ",
              np.count_nonzero(best_lr.feature_names_in_[(coefs[coefs != 0].dropna()).index]))

        fig, ax = plt.subplots()
        bar_coef = coefs.sort_values(by=0,ascending=False)[0:4].dropna()
        bar_coef.plot.bar(color="green")
        plt.xticks(range(0, len(best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=False)[0:4].dropna()).index])),
                   best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=False)[0:4].dropna()).index])
        plt.title(target_ab + ":\nStrongest positive coefficients for susceptibility prediction")
        plt.legend([])
        plt.tight_layout()
        plt.savefig(target_ab + targ_result + "-posfeatures.pdf", format="pdf", bbox_inches="tight")
        plt.show()

        fig, ax = plt.subplots()
        bar_coef = abs(coefs.sort_values(by=0,ascending=True)[0:4].dropna())
        bar_coef.plot.bar(color="green")
        plt.xticks(range(0, len(best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=True)[0:4].dropna()).index])),
                   best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=True)[0:4].dropna()).index])
        plt.title(target_ab + ":\nStrongest negative coefficients for susceptibility prediction")
        plt.legend([])
        plt.tight_layout()
        plt.savefig(target_ab + targ_result + "-negfeatures.pdf",format="pdf",bbox_inches="tight")
        plt.show()

    global c_value
    c_dict = logreg_cv.best_params_
    c_value = list(c_dict.values())[0]


#LR multinomial final run
def LR_multi_final(target,final_features,test=0.2,random=1,C_val=0.1,
             cw=None,model_type="auto",thr=0.5,
                   ab_of_interest="Ampicillin",av="micro",
                   targ_result="R"):

    global log_reg
    global fairness_model
    global probs
    probs = {}
    global class_report
    global model_dict
    model_dict = {}

    global X_test

    if len(target.unique()) > 2:

        for label in vari_list:

            X_train, X_test, y_train, y_test = train_test_split(
                final_features[vari_list[label]],
                target,
                test_size=test,
                random_state=random
            )

            # Train a new model and evaluate performance

            log_reg = LogisticRegression(
                max_iter=2000,
                multi_class=model_type,
                C=C_val,
                class_weight=cw,
                penalty='l1',
                solver='liblinear'
            )
            log_reg.fit(X_train, y_train)
            fairness_model = log_reg
            model_dict[label] = log_reg
            y_pred = (log_reg.predict_proba(X_test)[:, 1] >= thr).astype(int)
            y_pred_probs = log_reg.predict_proba(X_test)
            n_classes = len(np.unique(y_train))
            probs[label] = y_pred_probs

            # Generate ROC curve values: fpr, tpr, thresholds
            label_binarizer = LabelBinarizer().fit(y_train)
            y_onehot_test = label_binarizer.transform(y_test)
            y_onehot_test.shape  # (n_samples, n_classes)

           

            class_id = np.flatnonzero(label_binarizer.classes_ == label)[0]

            if label=="R":
                line_col="red"
            elif label=="I":
                line_col="darkorange"
            elif label=="S":
                line_col="green"
            elif label=="NT":
                line_col="blue"

            display = RocCurveDisplay.from_predictions(
                y_onehot_test[:, class_id],
                y_pred_probs[:, class_id],
                name=f"{label} vs the rest",
                color=line_col,
                plot_chance_level=True,
            )
            _ = display.ax_.set(
                xlabel="False Positive Rate",
                ylabel="True Positive Rate",
                title="One-vs-Rest ROC curve for "+ab_of_interest+":\n"+label+" vs rest",
            )

            plt.savefig(ab_of_interest + label + "-roc.pdf", format="pdf", bbox_inches="tight")
            plt.show()

        X_train, X_test, y_train, y_test = train_test_split(
            final_features[vari_list["R"]],
            target,
            test_size=test,
            random_state=random
        )

        # Train a new model and evaluate performance

        log_reg = LogisticRegression(
            max_iter=2000,
            multi_class=model_type,
            C=C_val,
            class_weight=cw,
            penalty='l1',
            solver='liblinear'
        )
        log_reg.fit(X_train, y_train)
        fairness_model = log_reg
        model_dict["R"] = log_reg
        y_pred = (log_reg.predict_proba(X_test)[:, 1] >= thr).astype(int)
        y_pred_probs = log_reg.predict_proba(X_test)
        n_classes = len(np.unique(y_train))
        classes = log_reg.classes_
        probs["R"] = y_pred_probs
        y_testpred = log_reg.predict(X_test)
        class_report = classification_report(y_test, y_testpred,output_dict=True)

        # Generate ROC curve values: fpr, tpr, thresholds
        label_binarizer = LabelBinarizer().fit(y_train)
        y_onehot_test = label_binarizer.transform(y_test)
        y_onehot_test.shape  # (n_samples, n_classes)

        display = RocCurveDisplay.from_predictions(
            y_onehot_test.ravel(),
            y_pred_probs.ravel(),
            name="micro-average OvR",
            color="purple",
            plot_chance_level=True,
        )
        _ = display.ax_.set(
            xlabel="False Positive Rate",
            ylabel="True Positive Rate",
            title=ab_of_interest + " Micro-averaged One-vs-Rest\nROC",
        )

        plt.savefig(ab_of_interest + "-avroc.pdf", format="pdf", bbox_inches="tight")
        plt.show()

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


        #plot all roc curves together

        fig, ax = plt.subplots(figsize=(6, 6))

        plt.plot(
            fpr["micro"],
            tpr["micro"],
            label=f"micro-average ROC curve (AUC = {roc_auc['micro']:.2f})",
            color="deeppink",
            linestyle=":",
            linewidth=4,
        )

        plt.plot(
            fpr["macro"],
            tpr["macro"],
            label=f"macro-average ROC curve (AUC = {roc_auc['macro']:.2f})",
            color="navy",
            linestyle=":",
            linewidth=4,
        )

        colors = cycle(["aqua", "darkorange", "cornflowerblue","limegreen"])
        for class_id, color in zip(range(n_classes), colors):
            RocCurveDisplay.from_predictions(
                y_onehot_test[:, class_id],
                y_pred_probs[:, class_id],
                name=f"ROC curve for {classes[class_id]}",
                color=color,
                ax=ax,
                plot_chance_level=(class_id == 2),
            )

        _ = ax.set(
            xlabel="False Positive Rate",
            ylabel="True Positive Rate",
            title=ab_of_interest+" ROCs for\nOne-vs-Rest multiclass",
        )

        plt.savefig(ab_of_interest + "-allroc.pdf", format="pdf", bbox_inches="tight")
        plt.show()

    else:


        X_train, X_test, y_train, y_test = train_test_split(
            final_features[vari_list],
            target,
            test_size=test,
            random_state=random
        )

        # Train a new model and evaluate performance

        log_reg = LogisticRegression(
            max_iter=2000,
            C=C_val,
            class_weight=cw,
            penalty='l1',
            solver='liblinear'
        )

        log_reg.fit(X_train, y_train)
        model_dict = log_reg
        fairness_model = log_reg
        y_pred = (log_reg.predict_proba(X_test)[:, 0] >= thr).astype(int)
        y_pred_probs = log_reg.predict_proba(X_test)[:,0]
        probs = y_pred_probs
        n_classes = len(np.unique(y_train))
        
        y_testpred = log_reg.predict(X_test)
        class_report = classification_report(y_test, y_testpred,output_dict=True)

        if log_reg.classes_[0] == "R":
            line_col = "red"
        elif log_reg.classes_[0] == "I":
            line_col = "darkorange"
        elif log_reg.classes_[0] == "S":
            line_col = "green"
        elif log_reg.classes_[0] == "NT":
            line_col = "blue"

        display = RocCurveDisplay.from_predictions(
            y_test,
            y_pred_probs,
            name=f"Binary prediction",
            color="purple",
            plot_chance_level=True,
            pos_label=log_reg.classes_[0]
        )
        _ = display.ax_.set(
            xlabel="False Positive Rate",
            ylabel="True Positive Rate",
            title="Binary ROC curve for " + ab_of_interest + " susceptibility"
        )

        plt.savefig(ab_of_interest + "-roc.pdf", format="pdf", bbox_inches="tight")
        plt.show()


def lr_hyp_tune_feats2(df,target,ht_features,test=0.2,random=1,folds=6,model_type="auto",cw=None,
                      scorer='f1_weighted',target_ab='antibiotic',targ_result=''):
    X_train, X_test, y_train, y_test = train_test_split(
        ht_features,
        target,
        test_size=test,
        random_state=random
    )
    log_reg = LogisticRegression(
        max_iter=2000,
        multi_class=model_type,
        class_weight=cw,
        penalty='l1',
        solver='liblinear'
    )
    param_grid = {"C": np.linspace(0.00001, 1, 10)}
    kf = KFold(n_splits=folds, random_state=random, shuffle=True)
    logreg_cv = GridSearchCV(log_reg, param_grid, cv=kf,scoring=scorer)
    logreg_cv.fit(X_train, y_train)
    print("Tuned logreg parameters: {}".format(logreg_cv.best_params_))
    print("Tuned logreg score: {}".format(logreg_cv.best_score_))
    n_scores = cross_val_score(log_reg, X_train, y_train, scoring=scorer, cv=kf, n_jobs=-1)
    print('Mean Accuracy: %.3f (%.3f)' % (stat.mean(n_scores), stat.stdev(n_scores)))
    best_lr = logreg_cv.best_estimator_
    global coefs
    coefs = pd.DataFrame(best_lr.coef_)
    coefs = coefs.transpose()
    if len(sorted(coefs)) > 1:
        coefs.columns = best_lr.classes_
    print("Total number of features:", len(ht_features.columns))
    global vari_list
    vari_list = {}

    if len(sorted(coefs)) > 1:
        for label in best_lr.classes_:
            vari_list[label] = best_lr.feature_names_in_[coefs[coefs[label] != 0].index]
            print("Nonzero features for " + label + ": ", vari_list[label])
            print("Number of selected features for " + label + ": ",
                  np.count_nonzero(best_lr.feature_names_in_[coefs[coefs[label] != 0].index]))

            if label == "R":
                line_col = "red"
            elif label == "I":
                line_col = "darkorange"
            elif label == "S":
                line_col = "green"
            elif label == "NT":
                line_col = "blue"


    else:
        vari_list = best_lr.feature_names_in_[(coefs[coefs != 0].dropna()).index]
        print("Nonzero features: ", vari_list)
        print("Number of selected features: ",
              np.count_nonzero(best_lr.feature_names_in_[(coefs[coefs != 0].dropna()).index]))

        fig, ax = plt.subplots()
        bar_coef = coefs.sort_values(by=0,ascending=False)[0:4].dropna()
        bar_coef.plot.bar(color="green")
        plt.xticks(range(0, len(best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=False)[0:4].dropna()).index])),
                   best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=False)[0:4].dropna()).index])
        plt.title(target_ab + ":\nStrongest positive coefficients for susceptibility prediction")
        plt.legend([])
        plt.tight_layout()
        plt.savefig(target_ab + targ_result + "-posfeatures.pdf", format="pdf", bbox_inches="tight")
        plt.show()

        fig, ax = plt.subplots()
        bar_coef = abs(coefs.sort_values(by=0,ascending=True)[0:4].dropna())
        bar_coef.plot.bar(color="green")
        plt.xticks(range(0, len(best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=True)[0:4].dropna()).index])),
                   best_lr.feature_names_in_[(coefs.sort_values(by=0,ascending=True)[0:4].dropna()).index])
        plt.title(target_ab + ":\nStrongest negative coefficients for susceptibility prediction")
        plt.legend([])
        plt.tight_layout()
        plt.savefig(target_ab + targ_result + "-negfeatures.pdf",format="pdf",bbox_inches="tight")
        plt.show()

    global c_value
    c_dict = logreg_cv.best_params_
    c_value = list(c_dict.values())[0]



#Model fairness assessment
def mod_fairness(protvarv,chosen_model,target,final_features,test=0.2,random=1,C_val=0.1,
             cw=None,model_type="auto",thr=0.5,
                   ab_of_interest="Ampicillin",av="micro",
                   targ_result="R",Bool=True):

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




#Model stability assessment

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


#iterated stability prediction

def iter_stab_pred(Antimicrobial_agent):

    global p1
    global p2
    global p3
    global p4
    global p5
    global p6
    global p7
    global p8
    global p9

    for i in list(range(1,101,1)):

        mod_stab_pred(target=urines5["Y"],
                      final_features=features,
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
        p2 = pd.concat([p2,pointtwo])

        pointthree = pd.DataFrame.from_dict(size_probs['0.86'])
        pointthree.columns = classes['0.86']
        pointthreeb = pd.DataFrame.from_dict(size_probs['actval0.86'])
        pointthree['actval'] = pointthreeb
        pointthree['AUC'] = size_aucs[0.86]
        pointthree['Antimicrobial'] =  Antimicrobial_agent
        pointthree['model'] = i
        p3 = pd.concat([p3,pointthree])

        pointfour = pd.DataFrame.from_dict(size_probs['0.88'])
        pointfour.columns = classes['0.88']
        pointfourb = pd.DataFrame.from_dict(size_probs['actval0.88'])
        pointfour['actval'] = pointfourb
        pointfour['AUC'] = size_aucs[0.88]
        pointfour['Antimicrobial'] =  Antimicrobial_agent
        pointfour['model'] = i
        p4 = pd.concat([p4,pointfour])

        pointfive = pd.DataFrame.from_dict(size_probs['0.9'])
        pointfive.columns = classes['0.9']
        pointfiveb = pd.DataFrame.from_dict(size_probs['actval0.9'])
        pointfive['actval'] = pointfiveb
        pointfive['AUC'] = size_aucs[0.9]
        pointfive['Antimicrobial'] =  Antimicrobial_agent
        pointfive['model'] = i
        p5 = pd.concat([p5,pointfive])

        pointsix = pd.DataFrame.from_dict(size_probs['0.92'])
        pointsix.columns = classes['0.92']
        pointsixb = pd.DataFrame.from_dict(size_probs['actval0.92'])
        pointsix['actval'] = pointsixb
        pointsix['AUC'] = size_aucs[0.92]
        pointsix['Antimicrobial'] =  Antimicrobial_agent
        pointsix['model'] = i
        p6 = pd.concat([p6,pointsix])

        pointseven = pd.DataFrame.from_dict(size_probs['0.94'])
        pointseven.columns = classes['0.94']
        pointsevenb = pd.DataFrame.from_dict(size_probs['actval0.94'])
        pointseven['actval'] = pointsevenb
        pointseven['AUC'] = size_aucs[0.94]
        pointseven['Antimicrobial'] =  Antimicrobial_agent
        pointseven['model'] = i
        p7 = pd.concat([p7,pointseven])

        pointeight = pd.DataFrame.from_dict(size_probs['0.96'])
        pointeight.columns = classes['0.96']
        pointeightb = pd.DataFrame.from_dict(size_probs['actval0.96'])
        pointeight['actval'] = pointeightb
        pointeight['AUC'] = size_aucs[0.96]
        pointeight['Antimicrobial'] =  Antimicrobial_agent
        pointeight['model'] = i
        p8 = pd.concat([p8,pointeight])

        pointnine = pd.DataFrame.from_dict(size_probs['0.98'])
        pointnine.columns = classes['0.98']
        pointnineb = pd.DataFrame.from_dict(size_probs['actval0.98'])
        pointnine['actval'] = pointnineb
        pointnine['AUC'] = size_aucs[0.98]
        pointnine['Antimicrobial'] =  Antimicrobial_agent
        pointnine['model'] = i
        p9 = pd.concat([p9,pointnine])

        print(i)

