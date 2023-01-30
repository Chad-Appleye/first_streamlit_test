import pickle
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
import os

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import mean_squared_error as MSE
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from sklearn.preprocessing import scale


class Training:
    """ 
    A class to train to machine learning models for CoPrimer predictions

    Attributes
    ----------
    data_path : str, path
        path to the screening data in an Excel workbook
    
    Methods
    ---------
    train_regression(output_dir)
        Trains the CoPrimer regression model for predicting Ct values

    train_classification(output_dir)
        Trains the primer dimer prediction model for predicting primer dimer formation in CoPrimer pairs

    train_all(output_dir)
        Trains both the regression and classification models


    """

    def __init__(self, data_path):
        # Train the classification and regression models on new data and verify the accuracy of the trained models
        #  param data: path to .XLSX containing training data in the correct format

        self.data_path = data_path
        self.seed = 42 #random number seed for reproducibility
        self.df_classification, self.df_regression = self._read_data()
        
    
    def _read_data(self):
        """Read the data and format it for the training of the regression and classification models
        Arguments:
            None
        Return:
            classification dataframe, regression dataframe"""

        df = pd.read_excel(self.data_path, sheet_name=0, engine='openpyxl')
        df.dropna(axis=1, how='all', inplace=True)
        df.dropna(axis=0, how='any', inplace=True)

        self.df = df

        #Drop the columns that don't contribute to the regression model
        df_regression = df.drop(columns=['PD Ct', 'Efficiency', '(f)OligoName', '(r)OligoName',
                      'Concatemer Formation'])

        #Drop the columns that don't contribute to the classification model
        df_classification = df.drop(columns=['Efficiency', '(f)OligoName', '(r)OligoName',
                      'Concatemer Formation', 'Ct '])

        #drop null values and get dummy variables
        df_regression = pd.get_dummies(df_regression)
        df_regression = df_regression.dropna()
        
        df_classification = pd.get_dummies(df_classification)
        df_classification = df_classification.dropna()

        #make a binary classification for primer dimer formation
        def binary(x):
            if x < 45:
                y = 1
            else:
                y = 0
            return y
        df_classification['PD'] = df_classification['PD Ct'].apply(binary)

        #drop PD CT column
        df_classification = df_classification.drop(columns='PD Ct')

        # exclude any Ct value < 15, as this value would not signify true amplification
        # of the target
        df_regression = df_regression[df['Ct '] > 15]

        return df_classification, df_regression

    
    def train_regression(self, output_dir):
        """Train only the regression model with the data specified

        Trains the regression model for prediction of Ct values of any given CoPrimer pair. The saved model and figure depicting the prediction evaluation will be saved in the output_dir
        Arguments:
            output_dir (str, path) -- directory to write the regression model 
        return: 
            None, pickled regression model file will be written in specified output_dir 
            prints evaluation of the trained model
        """

        #Seperate the data into features and target
        y = self.df_regression['Ct ']
        X = self.df_regression.drop(columns='Ct ')
        X_train, X_test, y_train, y_test = train_test_split(X, y,
                    test_size = .25, random_state=self.seed)
        #random forest regression
        
        rf = RandomForestRegressor(n_estimators=1577, oob_score=True, random_state=self.seed,
                           max_depth=50, max_features='auto', min_samples_leaf=1,
                           min_samples_split=2, n_jobs=-1)
    
        # fit the model to the training data
        rf.fit(X_train, y_train)

        #export trained model
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, 'regression_model.sav')
        pickle.dump(rf, open(output_file, 'wb'))

        y_pred = rf.predict(X_test)

        #compute the root mean squared error
        rmse_test = MSE(y_test, y_pred)**(1/2)
        SI = rmse_test/self.df_regression['Ct '].mean() 
        r2 = rf.score(X_test, y_test)

        print('Test set RMSE for Random Forest: {:.2f}'.format(rmse_test))
        print('Scatter Index for Random Forest: {}'.format(SI))
        print('Out of bag score: {}'.format(rf.oob_score_))#print out of bag score
        print('R^2 score: {}'.format(rf.score(X_test, y_test)))

        #sort the predicted and actual Ct values in pandas DF
        d = {'Actual Ct': y_test, 'Predicted Ct': y_pred}
        pred_df = pd.DataFrame(data=d)
        pred_df = pred_df.sort_values(by='Actual Ct')

        #plot the actual vs predicted Ct values
        self._scatter_fig(pred_df, r2, rmse_test, output_dir)


    def train_classification(self, output_dir):

        """ Train only the classification model

        Trains the binary classification model for primer dimer prediction of CoPrimer pairs. The saved model and figure depicting the prediction evaluation will be saved in the output_dir
        Arguments:
            output_dir (str, path) -- directory to write the regression model and evaluation
        return: 
            None, pickled classification model file will be written in specified output_dir 
            prints evaluation of the trained model
        """

        # Seperate the features and targets 
        y = self.df_classification['PD']
        X = self.df_classification.drop(columns='PD')

        # Split testing and training sets

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.25, random_state=self.seed)

        #classifier 
        rf_class = RandomForestClassifier(max_features='auto', min_samples_leaf=1, 
                            min_samples_split=2, n_estimators=700, 
                            oob_score=True)

        # Fit the model to the data
        rf_class.fit(X_train, y_train)

        # Export trained model
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, 'classification_model.sav')
        pickle.dump(rf_class, open(output_file, 'wb'))

        # Make predictions 
        y_pred = rf_class.predict(X_test)

        # Generate confusion matrix
        conf_mat = confusion_matrix(y_test, y_pred)
        accuracy = accuracy_score(y_test, y_pred)

        print('Accuracy score: {:.3f}'.format(accuracy))
        print('Out of bag score: {:.3f}'.format(rf_class.oob_score_))
        print('Confusion matrix: \n {}'.format(conf_mat))
        print('Classification Report:\n{}'.format(classification_report(y_test, y_pred)))

        #make dataframe of results split into actual vs predicted columns
        d = {'Actual': y_test, 'Predicted': y_pred}
        df_results = pd.DataFrame(d)

        self._bar_fig(df_results, output_dir, accuracy)


    def train_all(self, output_dir):
        """ Uses train_regression() and train_classification() methods to train both models
        Arguments:
            output_dir: (str, path) -- path to the output directory where models will be written
        Return:
            None. Models will be written and performance evaluations will be printed 
        """

        self.train_classification(output_dir=output_dir)
        self.train_regression(output_dir=output_dir)


    def _scatter_fig(self, df, r2, rmse, output_dir):
        """ Generate the scatter plot with plotly to evaluate the regression model """
            
        df.reset_index(drop=True, inplace=True)
        fig = go.Figure()

        fig.add_trace(go.Scatter(x=df.index, y=df['Actual Ct'], name='Actual', mode='markers', marker_color='blue'))
        fig.add_trace(go.Scatter(x=df.index, y=df['Predicted Ct'], name='Predicted', mode='markers', opacity=.7, marker_color='red'))
        fig.update_traces(marker=dict(size=5, line=dict(width=.3, color='DarkSlateGrey')))
        fig.update_layout(title_text='Predicted vs actual Ct Values', xaxis_title='Sample', yaxis_title='Ct Value')
        fig.add_annotation(x=100, y=50, text=' R^2 = {:.2}<br>RMSE: {:.2}'.format(r2, rmse), showarrow=False)

        fig.show()

        # write the file 
        fig_dir = os.path.join(output_dir, 'figures')
        os.makedirs(fig_dir, exist_ok=True)
        fig.write_html(os.path.join(fig_dir, 'regression_eval.html'))


    def _bar_fig(self, df, output_dir, accuracy):
        """ Generate the bar chart to evaluate the accuracy of the classification model  """

        df.reset_index(inplace=True, drop=True)
        df['Predicted_PD'] = df['Predicted'].apply(lambda x: 'Predicted' if x == 1 else 'Not Predicted')
        outcome = []
        for name, row in df.iterrows():
            if row['Actual'] == row['Predicted']:
                outcome.append('Correct')
            else:
                outcome.append('Incorrect')
        df['Outcome'] = outcome

        fig = px.histogram(df, x='Predicted_PD', color='Outcome', barmode='group')
        fig.add_annotation(x=0, y=len(df) - df['Predicted'].sum(), text='Accuracy score: {:.2}'.format(accuracy), showarrow=False)
        fig.update_layout(title_text='Primer Dimer Prediction Outcomes', xaxis_title='Predictions')
        fig.show()

        #write the file
        fig_dir = os.path.join(output_dir, 'figures')
        os.makedirs(fig_dir, exist_ok=True)
        fig.write_html(os.path.join(fig_dir, 'classification_eval.html'))


