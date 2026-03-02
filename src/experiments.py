import pandas as pd
from src.bridging import *

aptamers = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4,]
aptamers = [str(aptamer) for aptamer in aptamers]
df = pd.DataFrame(
    {
    'aptamer': aptamers,
    'gProcessedSignal': [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4,],
    'gIsFeatPopnOL': [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
    }
)
df1 = df.copy()

def aggregate(df1):
    df = df1.copy()
    for aptamer in df['aptamer'].unique():
        if aptamer[0].isdigit():
            df2 = df[df['aptamer']==aptamer]
            idx = df2.index
            relevant_df = df2[df2['gIsFeatPopnOL']!=1]
            mean = relevant_df['gProcessedSignal'].mean()
            df.loc[idx, 'gProcessedSignal'] = mean

    return df

df2 = aggregate(df1)
print()