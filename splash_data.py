import pandas as pd
import random

def find_outliers_IQR(df):
    q1=df.quantile(0.20)
    q3=df.quantile(0.80)
    IQR=q3-q1
    
    outliers = df[((df<(q1-1.5*IQR)) | (df>(q3+1.5*IQR)))]

    return outliers

high_speed_cam_df = pd.read_csv("high_speed_camera_1529.csv")
sticky_paper_df = pd.read_csv("sticky_paper_1529.csv")

# count splashes from 0
high_speed_cam_df['no'] -= 1
sticky_paper_df['no'] -= 1


# outliers handling

# print(find_outliers_IQR(high_speed_cam_c_df['v']).index.tolist())
# print(high_speed_cam_c_df['v'].max())
high_speed_cam_df = high_speed_cam_df.drop(find_outliers_IQR(high_speed_cam_df['v']).index.tolist())
# print(high_speed_cam_c_df['v'].max())



energy_sum_df = pd.DataFrame([high_speed_cam_df[high_speed_cam_df['no'] == i]['e'].sum() for i in range(16)], columns = ['energy_sum'])
# print(energy_sum_df)

energy_mean = energy_sum_df.mean()[0]
energy_std = energy_sum_df.std()[0]

print('Energy mean and std:')
print(energy_mean)
print(energy_std)

# cardinality of each splash
high_speed_cam_len = [len(high_speed_cam_df[high_speed_cam_df['no'] == i]) for i in range(16)]
sticky_paper_len = [len(sticky_paper_df[sticky_paper_df['no'] == i]) for i in range(48)]

high_speed_cam_len_sum = sum(high_speed_cam_len)
sp_hsc_card = pd.DataFrame([sum(random.sample(sticky_paper_len, 16))/sum(high_speed_cam_len) for i in range(1000)], columns = ['z'])

print(sp_hsc_card.mean(), sp_hsc_card.std())