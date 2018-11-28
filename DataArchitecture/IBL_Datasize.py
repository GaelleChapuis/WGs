import numpy as np
print('Training')
nlabs = np.array([10, 12]) # 10-12 experimental labs
n_mice_per_ses = np.array([4, 6]) # 4-6 mice at a time
ndays_train_per_year = round(365*5/7) # 261 days per year training
n_ses_per_year = ndays_train_per_year * nlabs * n_mice_per_ses
print(str(n_ses_per_year), ' training sessions per year')

# 2 Go/camera/hour of recording at 25 Hz, cropped, mpeg compression, 850/570 framesize, 2 cams
size_training_video_hourly_Gb = np.array([2,4])
ncams = 1
ses_duration_h = 1 # Session time 1 hour
size_training_videos_yearly_Tb = n_ses_per_year*ncams*size_training_video_hourly_Gb*ses_duration_h/1024
print(size_training_videos_yearly_Tb, 'Tb of training videos')

print('\n Recording')
mice_per_batch = np.array([1, 2]) # 1-2 mice recorded per batch
batch_cycle_length_days = 7*8 # 8 weeks from training to recording
n_mice_rec_per_year = np.round(365/batch_cycle_length_days*nlabs) # IBL mice recorded per year
rec_per_mouse = np.array([4, 6])  # recordings per mouse
n_rec_per_year = rec_per_mouse * n_mice_rec_per_year
print(n_rec_per_year , ' IBL recording sessions in a year')
rec_duration_h = 1.  # 1 hour per recording

ncams = np.array([3, 4]) # recording at 60 Hz, cropped, mpeg compression, 850/570 framesize
size_recording_video_hourly_Gb = np.array([4, 8])
size_recording_video_yearly_Tb = size_recording_video_hourly_Gb*n_rec_per_year*ncams*rec_duration_h/1024
print(size_recording_video_yearly_Tb, ' Tb of recording videos')

nprobes = 2
nchannels = 384
#Neuropixel 30kHz + 1kHz, 10bits, 384 channels, 2 probes
size_recording_neuropixel_hourly_Gb = 31000*nchannels*10/8*3600/1024/1024/1024*nprobes
size_recording_video_yearly_Tb = size_recording_neuropixel_hourly_Gb*n_rec_per_year*rec_duration_h/1024
print(size_recording_video_yearly_Tb , ' Tb of recording neuropixel data')


# Training
# [10440 18792]  training sessions per year
# [20.390625 36.703125] Tb of training videos
#  Recording
# [260. 468.]  IBL recording sessions in a year
# [4.0625 7.3125]  Tb of recording videos
# [25.33430234 45.6017442 ]  Tb of recording neuropixel data