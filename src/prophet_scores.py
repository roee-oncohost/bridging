import numpy as np


def rap_score_to_prophet_score(dev_rap_scores, subj_rap_score):
    prophet_score = 10*(dev_rap_scores<subj_rap_score).mean()
    return prophet_score

def rap_score_to_prophet_score_result(dev_rap_scores, subj_rap_score):
    prophet_score = 10*(dev_rap_scores<subj_rap_score).mean()
    prophet_result = (prophet_score >= 5).map({False: 'Negative', True: 'Positive'})
    return prophet_score, prophet_result

def get_prophet_score(df, model):
    df = df.copy()
    protein_columns = [col for col in df.columns if col[0].isdigit()]
    if df[protein_columns].max().max() > 50:
        df[protein_columns] = df[protein_columns].apply(np.log2)
    
    new_prediction_table = model.predict(df[protein_columns])
    rap_scores = new_prediction_table['y_pred_sp_scaled'].apply(lambda x: rap_score_to_prophet_score(model.prediction['y_pred_sp_scaled'], x))
    return rap_scores