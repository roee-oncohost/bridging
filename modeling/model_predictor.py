import json
import pickle
from os import path, sep, getcwd
from pandas.core.frame import DataFrame 
import numpy as np
import pandas as pd
try:
    from canopy.adat import Adat
except:
    from somadata.adat import Adat
from logger_manager import getLogger

logger = getLogger('RapResponsePredictor')

CONFIGURATION_FILE_NAME = 'response_prediction_configuration.json'


class ModelPredictor():
    def __init__(self):
        # load configuration file from script folder
        try:
            rap_predictor_working_dir = path.dirname(__file__) + sep
        except NameError:
            # Fallback to current working directory if __file__ is not defined
            rap_predictor_working_dir = getcwd() + sep
                
        # read configuration parameters
        with open(path.join(rap_predictor_working_dir, CONFIGURATION_FILE_NAME), 'r') as f:
            configuration_dict = json.load(f)['rap_response_predictor']
        rm        = configuration_dict['response_metrics'][0]
        model_file = configuration_dict[rm + '_initialization_params']['model_file']
        prophet_score_threshold = configuration_dict['prophet_score_threshold']
        
        # Load model
        self.model_file = path.join(rap_predictor_working_dir, model_file)
        self.prophet_score_threshold = prophet_score_threshold
        with open(self.model_file, 'rb') as f:
            self.model = pickle.load(f)
            
        self.pred_max = self.model.number_of_KS_proteins
        if type(self.model.prediction) == dict:
            self.dev_predictions = self.pred_max - pd.DataFrame.from_dict(self.model.prediction, orient='index')['y_pred_sp'].values
        elif type(self.model.prediction) == pd.DataFrame:
            self.dev_predictions = self.model.prediction['y_pred_sp'].values
        else:
            raise ValueError('Invalid type from model.prediction object')
        
        if hasattr(self.model, 'version'):
            self.reverse = False
        else:
            self.reverse = True

        ### End remove
        logger.info(f'{self.model_file} model loaded, dev SP AUC = {self.model.sp_auc:.6f}')
        
        
    def run_prediction(self, patient_data: DataFrame, sex: str = ''):
        ''' run prediction for a single patient, this method should be called by PRG '''
        if sex != '':
            valid_sex_vals = ['Male', 'Female', 'Other', 'Unknown', 'M', 'F', 'O', 'U']
            if len(patient_data) != 1:
                raise ValueError(f'Expecting a single-row adat - {len(patient_data)} does not match')
            if type(sex) != str or sex not in valid_sex_vals:
                raise ValueError(f'Expecting sex in str format (valid values: {", ".join(valid_sex_vals)})')
        
            if sex in ['Male', 'Female']:
                pass
            elif sex == 'M':
                sex = 'Male'
            elif sex == 'F':
                sex = 'Female'
            elif sex in ['Other', 'Unknown', 'O', 'U']:
                sex = 'Other'
            else:
                raise ValueError(f'"{sex}" is not a valid value for sex')
            
        patient_id = patient_data.index.get_level_values('SampleId')[0]
        logger.info(f'Running prediction for patient {patient_id}, sex = {sex}')
        if sex == 'Other':
            rv_male = self._run_prediction_multiple(patient_data, ['Male'])[patient_id]
            rv_fem  = self._run_prediction_multiple(patient_data, ['Female'])[patient_id]
            avg_pred = np.mean([rv_male['RAPScore'], rv_fem['RAPScore']])
            avg_prob = np.mean([rv_male['ResponseProbability'], rv_fem['ResponseProbability']])
            prpht_scr, prpht_res = self._rap_score_to_prophet_score(self.dev_predictions, avg_pred, self.prophet_score_threshold)
            rv = {'RAPScore': avg_pred,
                  'ResponseProbability': avg_prob,
                  'PROphetScore': prpht_scr,
                  'PROphetResult': prpht_res}
        elif sex == '':
            rv = self._run_prediction_multiple(patient_data, [])[patient_id]
        else:
            rv = self._run_prediction_multiple(patient_data, [sex])[patient_id]
       
        return rv
    
    
    def _rap_score_to_prophet_score(self, dev_rap_scores, subj_rap_score, prophet_score_threshold):
        prophet_score = (dev_rap_scores<subj_rap_score).mean()
        if prophet_score >= prophet_score_threshold:
            prophet_result = 'Positive'
        else:
            prophet_result = 'Negative'
        return prophet_score, prophet_result
    
    
    def _run_prediction_multiple(self, proteomics: Adat, sex: list = []):
        ''' run prediction for multiple patients '''
        if len(proteomics) != len(sex) and sex != []:
            raise ValueError('Proteomics Adat and list of patient sexes must be the same length')
    
        col_idx = [col for col in proteomics.columns.get_level_values('SeqId')] #[col + '_T0' for col in proteomics.columns.get_level_values('SeqId')]
        row_idx = proteomics.index.get_level_values('SampleId')
        patient_data = pd.DataFrame(np.log2(proteomics.values), index=row_idx, columns=col_idx)
        if sex != []:
            patient_data['Sex'] = pd.Series(sex, index=row_idx).map({'Male': 0, 'Female': 1})
    
        ### Remove if statement    
        if self.reverse:
            patient_data['Age'] = 60
            patient_data['PDL1'] = -1
            patient_data['Line'] = 0
        # Run model prediction for patient data
        pred = self.model.predict(patient_data)
        
        rv = dict()
        for subj in pred.index:
            ### Remove if statement
            if self.reverse:
                subj_pred = self.pred_max - pred.loc[subj, 'y_pred_sp']
                rv[subj] = {'RAPScore': float(subj_pred),
                            'ResponseProbability': float(1-pred.loc[subj, 'y_pred_sp_scaled']),
                            'PROphetScore': np.nan,
                            'PROphetResult': ''}
            else:
                subj_pred = pred.loc[subj, 'y_pred_sp']
                rv[subj] = {'RAPScore': float(subj_pred),
                            'ResponseProbability': float(pred.loc[subj, 'y_pred_sp_scaled']),
                            'PROphetScore': np.nan,
                            'PROphetResult': ''}
            
            # Calculate the response probability quantile relative to the dev set population
            rv[subj]['PROphetScore'], rv[subj]['PROphetResult'] = self._rap_score_to_prophet_score(self.dev_predictions, subj_pred, self.prophet_score_threshold)
            rv[subj]['PROphetScore'] = float(rv[subj]['PROphetScore'])
            
        return rv
    
if __name__ == '__main__':
    try:
        import canopy
    except:
        import somadata as canopy
    rcc_model_predictor = ModelPredictor()
    data_path = r"C:\Users\GilLoewenthal\Oncohost DX\Shares - Gil Loewenthal\code\sl_normalization\241027 Flag Pass check\014\OH2024_014.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat"
    # Deserialize the JSON string to a pandas DataFrame
    adat = canopy.read_adat(data_path)
    prediction = rcc_model_predictor.run_prediction(adat)