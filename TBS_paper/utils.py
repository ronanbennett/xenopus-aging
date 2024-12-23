import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import scipy
import sklearn.linear_model
import sklearn.model_selection
import sklearn.decomposition
from sklearn import set_config
set_config(enable_metadata_routing=True)
    
def make_data(train_meth_df, train_ages, test_meth_df, test_ages):
    X_train = train_meth_df.values.T
    y_train = np.array(train_ages)

    shared_names = [name for name in train_meth_df.index if name in set(test_meth_df.index)]
    X_train_adapted = train_meth_df.loc[shared_names].values.T
    y_train_adapted = y_train

    X_test = test_meth_df.loc[shared_names].values.T
    y_test = np.array(test_ages)

    print('X_train.shape', X_train.shape, 'y_train.shape', y_train.shape)
    print('X_train_adapted.shape', X_train_adapted.shape, 'y_train_adapted.shape', y_train_adapted.shape)
    print('X_test.shape', X_test.shape, 'y_test.shape', y_test.shape)
    return X_train, y_train, X_train_adapted, y_train_adapted, X_test, y_test


def clock(train_meth_df, train_ages, test_meth_df, test_ages, groups=None, transform=None, inv_transform=None, n_components=None, l1r_grid=None, L1R = 0.5, MAX_ITER = 8000):
    # X: np.array (n_samples, n_features)
    # y: np.array (n_samples, ) 
    X, y, X_train_adapted, y_train_adapted, X_test, y_test = make_data(train_meth_df, train_ages, test_meth_df, test_ages)
    assert(X.shape[0] == y.shape[0])

    is_predicting_transform = transform and inv_transform 
    if is_predicting_transform:
        y = transform(y)
        y_train_adapted = transform(y_train_adapted)
        y_test = transform(y_test)
    
    if groups is not None:
        outer_cv = sklearn.model_selection.LeaveOneGroupOut() #sklearn.model_selection.GroupKFold(n_splits=len(set(groups)))
        inner_cv = sklearn.model_selection.LeaveOneGroupOut()
    else:
        outer_cv = sklearn.model_selection.LeaveOneOut()  # cv = sklearn.model_selection.KFold(n_splits=192)
        inner_cv = 5
    
    # Apply PCA if n_components is specified
    if n_components:
        scaler = sklearn.preprocessing.StandardScaler()
        X = scaler.fit_transform(X)
        pca = sklearn.decomposition.PCA(n_components=n_components)
        X = pca.fit_transform(X)

    # # The following 2 lines are faster by ~5x, but I don't know of a way to run it where the .fit 
    # # method run on ElasticNetCV has a groups parameter that depends on the fold that is being run.
    # # This means that the inner cv must be something like 5 fold, which leaks tank information across internal 
    # # train/valid folds.
    # lm = sklearn.linear_model.ElasticNetCV(l1_ratio=L1R, max_iter=MAX_ITER, cv=inner_cv)
    # predictions = sklearn.model_selection.cross_val_predict(lm, X, y, groups=groups, cv=outer_cv, n_jobs=-1)

    predictions = np.zeros([y.shape[0]])
    for train_idx, val_idx in tqdm(outer_cv.split(X, groups=groups)):
        lm = sklearn.linear_model.ElasticNetCV(l1_ratio=L1R, max_iter=MAX_ITER, cv=inner_cv, n_jobs=-1)
        if groups is not None:
            lm.fit(X[train_idx], y[train_idx], groups=groups[train_idx])
        else:
            lm.fit(X[train_idx], y[train_idx])
        pred = lm.predict(X[val_idx])
        predictions[val_idx] = pred

    print(predictions)

    # test set
    if X_test is not None and y_test is not None:
        test_cv = inner_cv  # in order for the process that generated the CV plot above to be predictive of the process that generates the test predictions
        test_lm = sklearn.linear_model.ElasticNetCV(l1_ratio=L1R, max_iter=MAX_ITER, cv=test_cv, n_jobs=-1)
        if groups is not None:
            test_lm.fit(X_train_adapted, y_train_adapted, groups=groups)
        else:
            test_lm.fit(X_train_adapted, y_train_adapted)
        predictions_test = test_lm.predict(X_test)

    if is_predicting_transform:
        y = inv_transform(y)
        y_train_adapted = inv_transform(y_train_adapted)
        y_test = inv_transform(y_test)
        predictions = inv_transform(predictions)
        predictions_test = inv_transform(predictions_test)
    
    print(predictions_test.shape)
    print(predictions_test)

    plt.scatter(y, predictions, s=3, c='black')
    plt.xlabel('Actual Age (years)')
    plt.ylabel('Predicted Age (years)')
    plt.plot([min(y), max(y)], [min(y), max(y)], 'k--')     # TODO   #####

    mean_ae = sklearn.metrics.mean_absolute_error(y, predictions)
    median_ae = sklearn.metrics.median_absolute_error(y, predictions)
    rsq = sklearn.metrics.r2_score(y, predictions)
    me = sklearn.metrics.max_error(y, predictions)
    pearson_r = scipy.stats.pearsonr(predictions, y)[0]

    title_plot = f'''median_ae: {round(median_ae, 2)}
    R^2: {round(rsq, 2)}
    pearson r {pearson_r}
    mean_ae: {round(mean_ae, 2)}
    max_error: {round(me, 2)}
    alpha: 
    l1_ratio: {round(L1R, 5)}
    max_iter: {round(MAX_ITER, 5)}
    age_transform: {transform}
    leave_out_groups: {groups is not None}
    inner_cv {inner_cv}
    outer_cv {outer_cv}'''

    plt.title(title_plot)
    plt.show()

    
    predictions_train = predictions
    ######
    y = y_test
    predictions = predictions_test
    plt.scatter(y, predictions, s=3, c='black')
    plt.xlabel('Actual Age (years)')
    plt.ylabel('Predicted Age (years)')
    plt.plot([min(y), max(y)], [min(y), max(y)], 'k--')   # TODO   #####

    mean_ae = sklearn.metrics.mean_absolute_error(y, predictions)
    median_ae = sklearn.metrics.median_absolute_error(y, predictions)
    rsq = sklearn.metrics.r2_score(y, predictions)
    me = sklearn.metrics.max_error(y, predictions)
    pearson_r = scipy.stats.pearsonr(predictions, y)[0]

    title_plot = f'''median_ae: {round(median_ae, 2)}
    R^2: {round(rsq, 2)}
    pearson r {pearson_r}
    mean_ae: {round(mean_ae, 2)}
    max_error: {round(me, 2)}
    alpha: 
    l1_ratio: {round(L1R, 5)}
    max_iter: {round(MAX_ITER, 5)}
    age_transform: {transform}
    leave_out_groups: {groups is not None}
    inner_cv {inner_cv}
    outer_cv {outer_cv}
    test_cv {test_cv}'''

    plt.title(title_plot)
    plt.show()

    return predictions_train, predictions_test, lm, test_lm


def simple_clock(X, y, L1R=0.5, MAX_ITER=100000):
    # X: np.array (n_samples, n_features)
    # y: np.array (n_samples, ) 
    assert(X.shape[0] == y.shape[0])

    lm = sklearn.linear_model.ElasticNetCV(l1_ratio=L1R, max_iter=MAX_ITER, cv=5)
    cv = sklearn.model_selection.LeaveOneOut()
    predictions = sklearn.model_selection.cross_val_predict(lm, X, y, cv=cv, n_jobs=-1)

    # Plot
    plt.scatter(y, predictions, s=3, c='black')
    plt.xlabel('Actual Age')
    plt.ylabel('Predicted Age')
    plt.plot([min(y), max(y)], [min(y), max(y)], 'k--') 

    mean_ae = sklearn.metrics.mean_absolute_error(y, predictions)
    median_ae = sklearn.metrics.median_absolute_error(y, predictions)
    rsq = sklearn.metrics.r2_score(y, predictions)
    me = sklearn.metrics.max_error(y, predictions)
    pearson_r = scipy.stats.pearsonr(predictions, y)[0]

    title_plot = f'''median_ae: {round(median_ae, 2)}
    R^2: {round(rsq, 2)}
    pearson r {pearson_r}
    pearson r, squared {pearson_r ** 2}
    mean_ae: {round(mean_ae, 2)}
    max_error: {round(me, 2)}
    l1_ratio: {round(L1R, 5)}
    max_iter: {round(MAX_ITER, 5)}'''

    plt.title(title_plot)
    plt.show()

    return predictions, cv, lm