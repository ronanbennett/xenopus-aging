import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, RegressorMixin

class XTropicalisClock(BaseEstimator, RegressorMixin):
    def __init__(self, coef_path):
        # Initializes  clock with coefficients from TSV file.
        self.coef_path = coef_path
        self._load_coefficients()

    def _load_coefficients(self):
        self.coef_df = pd.read_csv(self.coef_path, sep='\t')
        self.intercept_ = self.coef_df[self.coef_df['Position'] == 'Intercept']['Coefficient'].values[0]
        self.coef_df = self.coef_df[self.coef_df['Position'] != 'Intercept']
        self.feature_names_ = self.coef_df.apply(lambda row: f"{row['Chr']}:{row['Position']}", axis=1).values
        self.coef_ = self.coef_df['Coefficient'].values

    def fit(self, X, y=None):
        return self

    def predict(self, X):
        # X: pd.DataFrame, Rows are loci in 'Chr:Position', Columns are samples.
        # returns np.ndarray Predicted ages for each sample.
        common_features = list(set(self.feature_names_).intersection(set(X.index)))
        if not common_features:
            raise ValueError("No overlapping features between methylation data and model coefficients.")
        # Align the features according to the model's expected order
        X_ordered = X.loc[self.feature_names_].values.T
        y_pred = np.dot(X_ordered, self.coef_) + self.intercept_
        return y_pred
