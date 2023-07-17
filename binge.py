## functions for binge eating project

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.ticker import MaxNLocator
from matplotlib import gridspec
import seaborn as sns
import os
import glob
import datetime
import numpy as np
import pandas as pd
import math
from random import randrange
from random import randint
from tqdm import tqdm
import dabest
import sys
from BaselineRemoval import BaselineRemoval
import pynapple as nap
from pathlib import Path


################## IO ##########################
################################################
################## phenosys ####################

## Def phensys converter for lick rate
def timestampconvert(x):
  stamp = datetime.timedelta(days = x)
  result = datetime.datetime(1899,12,30,0,0) + stamp
  # print(result.strftime('%Y-%m-%d %H:%M:%S.%f'))
  return result

def datetime_convert_phenosys(csv_path):

    df_pheno = pd.read_csv(csv_path,sep=',')
    time_list = df_pheno['DateTime']
    result = []
    for x in time_list:
        timestampconvert(float(x))
        result.append(timestampconvert(x))
    new_time_stamps = []
    for i in result:
        new_time_stamps.append((i-min(result)).total_seconds())
    df_pheno = df_pheno.fillna(0)
    df_pheno['DateTime'] = new_time_stamps
    # df_pheno
    return df_pheno

## Function for getting L1 and L2 lick events
def lick_event_calculate(csv_path):
    """
    This function will calcuate the lick sensor data from the phenosys csv files and return 3 lists of timestamps from sensor L1, L2 and the timestamps of both channels. It will also do a quick plotting for the lick sensor data for raster plots and the density plot of overall lick events.

    """
    df = datetime_convert_phenosys(csv_path)
    
    # Get L1 timestamps and MsgValue1
    L1_time = df[df['unitLabel']=='L1']['DateTime'].values
    L1_value = df[df['unitLabel']=='L1']['MsgValue1'].values

    L1_timestamps_new = []
    for idx,value in enumerate(L1_value):
        if len(value.split(','))>1:
            for j in value.split(',')[1:]:
                individual_value = int(j.split('-')[0])*0.001 # convert ms to s
                if idx == 0:
                    L1_timestamps_new.append(L1_time[idx]+individual_value)
                else:
                    L1_timestamps_new.append(L1_timestamps_new[-1]+individual_value)
        else:
            L1_timestamps_new.append(L1_time[idx])
    L1_timestamps_new = np.array(L1_timestamps_new)

    ## The same for L2
    L2_time = df[df['unitLabel']=='L2']['DateTime'].values
    L2_value = df[df['unitLabel']=='L2']['MsgValue1'].values

    L2_timestamps_new = []
    for idx,value in enumerate(L2_value):
        if len(value.split(','))>1:
            for j in value.split(',')[1:]:
                individual_value = int(j.split('-')[0])*0.001
                if idx == 0:
                    L2_timestamps_new.append(L2_time[idx]+individual_value)
                else:
                    L2_timestamps_new.append(L2_timestamps_new[-1]+individual_value)
        else:
            L2_timestamps_new.append(L2_time[idx])
    L2_timestamps_new = np.array(L2_timestamps_new)

    all_lick_events = np.array(sorted([*L1_timestamps_new,*L2_timestamps_new]))

    return L1_timestamps_new, L2_timestamps_new, all_lick_events


## Function for getting Pump events
def Pump_event_calculate(csv_path):
    """
    This function will calcuate the lick sensor data from the phenosys csv files and return 3 lists of timestamps from GPIOs P1 and P1C and the timestamps of both channels.

    """
    df = datetime_convert_phenosys(csv_path)
    
    # Get L1 timestamps and MsgValue1
    if 'P1' in df['unitLabel'].unique():
      P1_time = df[df['unitLabel']=='P1']['DateTime'].values
    else:
      P1_time = df[df['unitLabel']=='P1D']['DateTime'].values ## I later switched to P1D pump since Pump1 is dead

    P1C_time = df[df['unitLabel']=='P1C']['DateTime'].values
    
    P2_time = df[df['unitLabel']=='P2']['DateTime'].values
    P2A_time = df[df['unitLabel']=='P2A']['DateTime'].values

    # all_pump_events = P1_time + P1C_time + P2_time + P2A_time
    pump_all = np.concatenate([P1_time,P1C_time,P2_time,P2A_time],axis=0)
    pump_all = sorted(pump_all)
    pump_all = np.array(pump_all)

    # fig, ax0 = plt.subplots(nrows=1,ncols=1,sharex=True,figsize=[8,2])
    # ax0.eventplot([P1_time,P1C_time,P2_time,P2A_time],lw=0.5,linelengths=0.8,color=['C0','C1','C3','C4'])
    # legend=ax0.legend(['P1','P1C','P2','P2A'],bbox_to_anchor=(0., 1.1, 1., 1.1), loc=3, ncol=4, mode="expand", 
    #             borderaxespad=0.,frameon=False,title='GPIO events: '+csv_path.split('/')[-1],fontsize=8)
    
    # legend.get_title().set_fontsize('9')
    
    # ax0.set_yticks([])
    # plt.xlabel('Second',fontsize=8)
    # plt.xticks(fontsize=8);plt.yticks(fontsize=8)
    # plt.xlim([0,1600])
    # plt.tight_layout();plt.show()
    
    return P1_time, P1C_time, P2_time, P2A_time, pump_all

def add_missing_licks(Pump_all_no_init,L1_timestamps):
    L1_timestamps = np.concatenate([np.zeros(1),L1_timestamps]) ## add a 0 at the beginning, so we don't need a conditional statement

    # diff_concat = []
    for pump_event in Pump_all_no_init:
        # print(pump_event)
        # if (L1_timestamps<pump_event).all() ==0:
        #    diff_ = 100 # set to a large value that passes the criteria
        # else:
        nearest_event = find_nearest(L1_timestamps[L1_timestamps<pump_event], pump_event) # only check the licks preceding the pump event
        # print(nearest_event,pump_event)
        diff_ = pump_event - nearest_event
        # print(diff_)
        # diff_concat.append(diff_)
        # if diff_<0: # means the nearest lick is after pump activation
        #     L1_timestamps = np.append(L1_timestamps,pump_event-0.02) ## mimicking the missing lick events
        if diff_>0.03: # measn neatest lick is way before pump activation
            # print('yes')
            L1_timestamps = np.append(L1_timestamps,pump_event-0.02) ## mimicking the missing lick events
        # else:
        #    print('lick event is presented, skip')
    L1_timestamps = np.array(sorted(L1_timestamps)[1:]) # remove the first 0
    return L1_timestamps

def detect_feeding_bout(Pump1,interval=3):
    # if pump interval is less than 2 seconds?
    Pump1 = np.pad(Pump1, (1, 0), 'constant', constant_values=(0,0))
    diff = Pump1[1:]-Pump1[:-1]
    gap = diff>interval # binary data
    idx = np.where(gap == 1) # idx == 1
    start_concat, end_concat = [],[]
    for idx_ in idx:
        start_concat.append(Pump1[idx_+1])
        end_concat.append(Pump1[idx_])
    start_concat = start_concat[0][:-1] # remove empty starting
    end_concat = end_concat[0][1:] # remove empty ending
    
    return start_concat, end_concat

def baseline_finish_time(csv_path):
    """
    This function will calcuate the lick sensor data from the phenosys csv files and return 3 lists of timestamps from sensor L1, L2 and the timestamps of both channels. It will also do a quick plotting for the lick sensor data for raster plots and the density plot of overall lick events.

    """
    df = datetime_convert_phenosys(csv_path)
    
    # Get L1 timestamps and MsgValue1
    baseline_finish = df[df['SystemMsg']=='baselinelickfinished']['DateTime'].values
    return baseline_finish[0]

# def remove_pheno_init(Pump1,Pump2):
#   """
#   this is for the inscopix imaging protocol
#   so there is 5 min baseline, and then 2 min/1 min feeding bout/interval

#   """
#   pump_all = [*Pump1,*Pump2]
#   pump_all = sorted(pump_all)
#   pump_all = np.array(pump_all)
#   x_shift = pump_all[-1]
#   x_bar = [x*180+300 for x in range(18)]
#   # x_bar = [0] + x_bar
#   x_bar = np.array(x_bar)
#   x_bar = x_bar[x_bar<x_shift]
#   for x in x_bar:
#     if len(Pump2) == 0: ## again to prevent the 2nd pump is empty
#       nearest_num = [find_nearest(Pump1,x)]
#       if find_nearest(nearest_num,x) in Pump1:
#         Pump1 = np.delete(Pump1, np.where(Pump1 == find_nearest(nearest_num,x)))
#     else:
#       nearest_num = [find_nearest(Pump1,x),find_nearest(Pump2,x)]
#       if find_nearest(nearest_num,x) in Pump1:
#         Pump1 = np.delete(Pump1, np.where(Pump1 == find_nearest(nearest_num,x)))
#       elif find_nearest(nearest_num,x) in Pump2:    
#         Pump2 = np.delete(Pump2, np.where(Pump2 == find_nearest(nearest_num,x)))
#   Pump_all_no_init = [*Pump1,*Pump2]
#   Pump_all_no_init = np.array(sorted(Pump_all_no_init))
#   return Pump1,Pump2,Pump_all_no_init

def remove_init(Pump1,Pump2):
  pump_all = [*Pump1,*Pump2]
  pump_all = sorted(pump_all)
  pump_all = np.array(pump_all)
  x_shift = pump_all[-1]
  x_bar = [x*180+300 for x in range(18)]
  x_bar = [0] + x_bar
  x_bar = np.array(x_bar)
  x_bar = x_bar[x_bar<x_shift]
  for x in x_bar:
    if len(Pump2) == 0: # again to prevent the 2nd pump is empty
      nearest_num = [find_nearest(Pump1,x)]
      if find_nearest(nearest_num,x) in Pump1:
        Pump1 = np.delete(Pump1, np.where(Pump1 == find_nearest(nearest_num,x)))
    else:
      if len(Pump1) ==0:
          if len(Pump2) == 1:
              Pump2 = []
          else:
              nearest_num = find_nearest(Pump2,x)
          if find_nearest(nearest_num,x) in Pump2: 
            Pump2 = np.delete(Pump2, np.where(Pump2 == find_nearest(nearest_num,x)))
      else:
          nearest_num = [find_nearest(Pump1,x),find_nearest(Pump2,x)]
          if find_nearest(nearest_num,x) in Pump1:
              Pump1 = np.delete(Pump1, np.where(Pump1 == find_nearest(nearest_num,x)))
          elif find_nearest(nearest_num,x) in Pump2:    
              Pump2 = np.delete(Pump2, np.where(Pump2 == find_nearest(nearest_num,x)))
    
  Pump_all_no_init = [*Pump1,*Pump2]

  if len(Pump_all_no_init) >0:
      Pump_all_no_init = np.array(sorted(Pump_all_no_init))
  return Pump1,Pump2,Pump_all_no_init



def load_phenosys(file_path):
    """
    Load following measures from a given phenosys csv file
    L1,L2:          (np.array) Lick senseor timestamps
    P1,P2,PAll:     (np.array) Pump events for milk & water & all together
    """
    ## load data
    L1, L2, all_lick_events  = lick_event_calculate(csv_path=file_path)
    P1_time, P1C_time, P2_time, P2A_time, pump_all = Pump_event_calculate(csv_path=file_path)

    L1 = select_lick_sensor(L1, L2)
    P1, P2  = select_pump(P1_time, P1C_time, P2_time, P2A_time)
    ## remove 0s in P1, P2
    P1 = [i for i in P1 if i >=300]
    P2 = [i for i in P2 if i >=300]
    # P1,P2,PAll = remove_pheno_init(Pump1,Pump2)
    if '21.01.06' in file_path and 'DSC012849' in file_path:
      P1,P2,PAll = remove_GPIO_init_spec(P1,P2)
    else:
      P1,P2,PAll = remove_GPIO_init(P1,P2)         

    L1= add_missing_licks(PAll,L1) ## notice I used Pump1_new, not Pump_all_no_init, since Pump2 data is not actually related to lick events
    # plt.eventplot(P1)
    return L1,P1,P2,PAll

def load_phenosys_nofilter(file_path):
    """
    Load following measures from a given phenosys csv file
    L1,L2:          (np.array) Lick senseor timestamps
    P1,P2,PAll:     (np.array) Pump events for milk & water & all together
    """
    ## load data
    L1, L2, all_lick_events  = lick_event_calculate(csv_path=file_path)
    P1_time, P1C_time, P2_time, P2A_time, PAll = Pump_event_calculate(csv_path=file_path)

    L1 = select_lick_sensor(L1, L2)
    P1, P2  = select_pump(P1_time, P1C_time, P2_time, P2A_time)
    ## remove 0 from pump event lists
    P1 = [i for i in P1 if i >=300]
    P2 = [i for i in P2 if i >=300]
    PAll = sorted([*P1,*P2])
    PAll = [i for i in PAll if i != 0]
    # P1,P2,PAll = remove_GPIO_init(P1,P2)
    # print(len(P1),len(P2),len(L1))
    # print(PAll)

    L1= add_missing_licks(PAll,L1) ## notice I used Pump1_new, not Pump_all_no_init, since Pump2 data is not actually related to lick events
    # plt.eventplot(P1)
    return L1,P1,P2,PAll

def pheno_eventplot(P1,P2,PAll,P1_,P2_,PAll_,L1,mouse_id,date,time_difference):   
    L1 = nap.Ts(L1)
    lick_max = np.max(L1.count(1).values)
    # print(lick_max)
    if len(PAll_) == 0:
      x_shift = 1550
    else:
      x_shift = PAll_[-1]
    if x_shift< 1550:
       x_shift = 1550 # force output to be 8 panels for consistency
    x_bar = [x*180+300 for x in range(18)]
    x_bar = [0] + x_bar
    x_bar = np.array(x_bar)
    x_bar = x_bar[x_bar<x_shift]
    n=0
    plt.figure(figsize=[12,4])
    plt.suptitle(f'{mouse_id}_{date}')
    for x in x_bar:
      ax = plt.subplot(2, (len(x_bar)+1)//2, n+1)
      ax.text(x=x,y=5,s=lick_max,ha='left',va='top')
      ax.plot(L1.count(1).index,(L1.count(1).values)/lick_max+4,alpha=0.75,c='grey')
      ax.eventplot([P1,P2,P1_,P2_],linewidths=1,colors=['C0','C0','C1','C1'],linelengths=0.8)
      ax.plot([P1[0],P1[0]-time_difference],[2,2],lw=2)
      ax.set_xlim([x-10,x+130])
      ax.set_ylim([-1,5])
      sns.despine(ax=ax)
      ax.set_yticks([])
      n+=1
    plt.tight_layout()
    plt.show()

def pheno_lick_eventplot(L1,L1_,PAll,mouse_id,date):   
    L1 = nap.Ts(L1)
    lick_max = np.max(L1.count(1).values)
    L1_ = nap.Ts(L1_)

    # print(lick_max)
    x_shift = L1.count(1).index[-1]
    if x_shift< 1550:
       x_shift = 1550 # force output to be 8 panels for consistency
    x_bar = [x*180+300 for x in range(18)]
    x_bar = [0] + x_bar
    x_bar = np.array(x_bar)
    x_bar = x_bar[x_bar<x_shift]
    n=0
    plt.figure(figsize=[12,4])
    plt.suptitle(f'{mouse_id}_{date}')
    for x in x_bar:
      ax = plt.subplot(2, (len(x_bar)+1)//2, n+1)
      ax.text(x=x,y=3,s=lick_max,ha='left',va='top')
      ax.plot(L1.count(1).index,(L1.count(1).values)/lick_max+1,alpha=1,c='C0')
      ax.plot(L1_.count(1).index,(L1_.count(1).values)/lick_max+2,alpha=1,c='C1')
      ax.eventplot(PAll,linelengths=0.8,linewidths=1,lineoffsets=-0.25)
      ax.set_xlim([x-10,x+130])
      # ax.set_ylim([-1,5])
      sns.despine(ax=ax)
      ax.set_yticks([])
      n+=1
    plt.tight_layout()
    plt.show()

## Filter function for np array
def filter_mask_large(arr, k):
    return arr[arr < k]
def filter_mask_small(arr, j):
    return arr[arr > j]

def filter_mask_range(arr, min, max):
    arr = filter_mask_large(arr, max)
    arr = filter_mask_small(arr, min)
    return arr

## find nearst value in np array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def select_pump(P1,P1C,P2,P2A):
  if (len(P1) + len(P1C)) > (len(P2) + len(P2A)):
    # print('Select P1/P1C')
    Pump1 = P1
    Pump2 = P1C
  else:
    # print('Select P2/2A')
    Pump1 = P2
    Pump2 = P2A
  return Pump1, Pump2

def select_lick_sensor(L1_timestamps, L2_timestamps):
  if len(L1_timestamps)>len(L2_timestamps):
    licks = L1_timestamps
  else:
    licks = L2_timestamps
  return licks

def select_opto_pump_led(P1,P1C,P2,P2A):
  if len(P1) > len(P2):
    Pump1 = P1
    Pump2 = P2
  else:
    Pump1 = P2
    Pump2 = P1
  return Pump1, Pump2

def remove_init(Pump1,Pump2):
  pump_all = [*Pump1,*Pump2]
  pump_all = sorted(pump_all)
  pump_all = np.array(pump_all)
  x_shift = pump_all[-1]
  x_bar = [x*180+300 for x in range(18)]
  x_bar = [0] + x_bar
  x_bar = np.array(x_bar)
  x_bar = x_bar[x_bar<x_shift]
  for x in x_bar:
    if len(Pump2) == 0: ## to prevent the 2nd pump is empty
      if len(Pump1) !=0:
        nearest_num = [find_nearest(Pump1,x)]
        if find_nearest(nearest_num,x) in Pump1:
          Pump1 = np.delete(Pump1, np.where(Pump1 == find_nearest(nearest_num,x)))
    else:
      if len(Pump1) != 0:
        nearest_num = [find_nearest(Pump1,x),find_nearest(Pump2,x)]
        if find_nearest(nearest_num,x) in Pump1:
          Pump1 = np.delete(Pump1, np.where(Pump1 == find_nearest(nearest_num,x)))
        elif find_nearest(nearest_num,x) in Pump2:    
          Pump2 = np.delete(Pump2, np.where(Pump2 == find_nearest(nearest_num,x)))
  Pump_all_no_init = [*Pump1,*Pump2]
  if len(Pump_all_no_init) == 0:
    Pump_all_no_init = np.array([])
  else:    
    Pump_all_no_init = np.array(sorted(Pump_all_no_init))
  return Pump1,Pump2,Pump_all_no_init

def remove_pheno_init(Pump1,Pump2):
  pump_all = [*Pump1,*Pump2]
  pump_all = sorted(pump_all)
  pump_all = np.array(pump_all)

  x = 0
  ## pump1 remove init
  nearest_num = [find_nearest(Pump1,x)]
  if find_nearest(nearest_num,x) in Pump1:
    Pump1 = np.delete(Pump1, np.where(Pump1 == find_nearest(nearest_num,x)))
  ## pump2 remove init
  nearest_num = [find_nearest(Pump2,x)]
  if find_nearest(nearest_num,x) in Pump2:
    Pump2 = np.delete(Pump2, np.where(Pump2 == find_nearest(nearest_num,x)))

  Pump_all_no_init = [*Pump1,*Pump2]
  if len(Pump_all_no_init) == 0:
    Pump_all_no_init = np.array([])
  else:    
    Pump_all_no_init = np.array(sorted(Pump_all_no_init))
  return Pump1,Pump2,Pump_all_no_init




def remove_GPIO_init(Pump1,Pump2):
  pump_all = [*Pump1,*Pump2]
  pump_all = sorted(pump_all)
  ## remove all pump event in 300 sec, since there shouldn't be any
  pump_all = [i for i in pump_all if i >= 300]
  pump_all = np.array(pump_all)
  x_shift_end = pump_all[-1]
  x_shift_start = pump_all[0] -120 # force the x_bar to start at the beginning of this round
  x_bar = [x*180+300 for x in range(18)]
  # x_bar = [0] + x_bar
  x_bar = np.array(x_bar)
  x_bar = x_bar[x_bar<x_shift_end]
  x_bar = x_bar[x_bar>x_shift_start]
  Pump1 = np.array(Pump1)
  Pump2 = np.array(Pump2)
  for x in x_bar:
    if find_nearest(pump_all[pump_all>x],x)-x>3: # this means there is no initial pulse within 3 sec at the round starting time
       pass
    else:
      if len(Pump2[Pump2>x]) == 0: ## prevent the 2nd pump is empty
        nearest_num = [find_nearest(Pump1[Pump1>x],x)]
        if find_nearest(nearest_num,x) in Pump1:
          Pump1 = np.delete(Pump1, np.where(Pump1 == find_nearest(nearest_num,x)))
          # print(f'{x}, del {find_nearest(nearest_num,x)} from Pump1')
      else:
        if len(Pump1[Pump1>x]) == 0:
          nearest_num = [find_nearest(Pump2[Pump2>x],x)]
          Pump2 = np.delete(Pump2, np.where(Pump2 == find_nearest(nearest_num,x)))
          # print(f'{x}, del {find_nearest(nearest_num,x)} from Pump2')
        else:
          nearest_num = [find_nearest(Pump1[Pump1>x],x),find_nearest(Pump2[Pump2>x],x)]
          if find_nearest(nearest_num,x) in Pump1:
            Pump1 = np.delete(Pump1, np.where(Pump1 == find_nearest(nearest_num,x)))
            # print(f'{x}, del {find_nearest(nearest_num,x)} from Pump1')
          elif find_nearest(nearest_num,x) in Pump2:    
            Pump2 = np.delete(Pump2, np.where(Pump2 == find_nearest(nearest_num,x)))
            # print(f'{x}, del {find_nearest(nearest_num,x)} from Pump2')
  Pump_all_no_init = [*Pump1,*Pump2]
  Pump_all_no_init = np.array(sorted(Pump_all_no_init))
  return Pump1,Pump2,Pump_all_no_init

def remove_isx_GPIO_init(Pump1,Pump2):
  pump_all = [*Pump1,*Pump2]
  pump_all = sorted(pump_all)
  ## remove all pump event in 300 sec, since there shouldn't be any
  pump_all = [i for i in pump_all if i >= 300]
  pump_all = np.array(pump_all)
  x_shift_end = pump_all[-1]
  x_shift_start = pump_all[0] -120 # force the x_bar to start at the beginning of this round
  x_bar = [x*180+300 for x in range(18)]
  # x_bar = [0] + x_bar
  x_bar = np.array(x_bar)
  x_bar = x_bar[x_bar<x_shift_end]
  x_bar = x_bar[x_bar>x_shift_start]
  Pump1 = np.array(Pump1)
  Pump2 = np.array(Pump2)
  for x in x_bar:
    if find_nearest(pump_all[pump_all>x],x)-x>25: # this means there is no initial pulse within 3 sec at the round starting time
       pass
    else:
      if len(Pump2[Pump2>x]) == 0: ## prevent the 2nd pump is empty
        if len(Pump1[Pump1>x]) == 0:
          #  print('no event in this range, skip')
          pass
        else:
          nearest_num = [find_nearest(Pump1[Pump1>x],x)]
          if find_nearest(nearest_num,x) in Pump1:
            Pump1 = np.delete(Pump1, np.where(Pump1 == find_nearest(nearest_num,x)))
            # print(f'{x}, del {find_nearest(nearest_num,x)} from Pump1 isx')
      else:
        if len(Pump1[Pump1>x]) == 0:
          nearest_num = [find_nearest(Pump2[Pump2>x],x)]
          Pump2 = np.delete(Pump2, np.where(Pump2 == find_nearest(nearest_num,x)))
          # print(f'{x}, del {find_nearest(nearest_num,x)} from Pump2 isx')
        else:
          nearest_num = [find_nearest(Pump1[Pump1>x],x),find_nearest(Pump2[Pump2>x],x)]
          if find_nearest(nearest_num,x) in Pump1:
            Pump1 = np.delete(Pump1, np.where(Pump1 == find_nearest(nearest_num,x)))
            # print(f'{x}, del {find_nearest(nearest_num,x)} from Pump1 isx')
          elif find_nearest(nearest_num,x) in Pump2:    
            Pump2 = np.delete(Pump2, np.where(Pump2 == find_nearest(nearest_num,x)))
            # print(f'{x}, del {find_nearest(nearest_num,x)} from Pump2 isx')
  Pump_all_no_init = [*Pump1,*Pump2]
  Pump_all_no_init = np.array(sorted(Pump_all_no_init))
  return Pump1,Pump2,Pump_all_no_init


def remove_GPIO_init_spec(Pump1,Pump2):
  pump_all = [*Pump1,*Pump2]
  pump_all = sorted(pump_all)
  ## remove all pump event in 300 sec, since there shouldn't be any
  pump_all = [i for i in pump_all if i >= 300]
  pump_all = np.array(pump_all)
  x_shift_end = pump_all[-1]
  x_shift_start = pump_all[0] -120 # force the x_bar to start at the beginning of this round
  x_bar = [x*180+300 for x in range(18)]
  # x_bar = [0] + x_bar
  x_bar = np.array(x_bar)
  x_bar = x_bar[x_bar<x_shift_end]
  x_bar = x_bar[x_bar>x_shift_start]
  x_bar = x_bar[0:-1]
  x_bar = np.array(list(x_bar) +[2032])
  Pump1 = np.array(Pump1)
  Pump2 = np.array(Pump2)
  for x in x_bar:
    if find_nearest(pump_all[pump_all>x],x)-x>3: # this means there is no initial pulse within 3 sec at the round starting time
       pass
    else:
      if len(Pump2[Pump2>x]) == 0: ## prevent the 2nd pump is empty
        nearest_num = [find_nearest(Pump1[Pump1>x],x)]
        if find_nearest(nearest_num,x) in Pump1:
          Pump1 = np.delete(Pump1, np.where(Pump1 == find_nearest(nearest_num,x)))
          # print(f'{x}, del {find_nearest(nearest_num,x)} from Pump1')
      else:
        if len(Pump1[Pump1>x]) == 0:
          nearest_num = [find_nearest(Pump2[Pump2>x],x)]
          Pump2 = np.delete(Pump2, np.where(Pump2 == find_nearest(nearest_num,x)))
          # print(f'{x}, del {find_nearest(nearest_num,x)} from Pump2')
        else:
          nearest_num = [find_nearest(Pump1[Pump1>x],x),find_nearest(Pump2[Pump2>x],x)]
          if find_nearest(nearest_num,x) in Pump1:
            Pump1 = np.delete(Pump1, np.where(Pump1 == find_nearest(nearest_num,x)))
            # print(f'{x}, del {find_nearest(nearest_num,x)} from Pump1')
          elif find_nearest(nearest_num,x) in Pump2:    
            Pump2 = np.delete(Pump2, np.where(Pump2 == find_nearest(nearest_num,x)))
            # print(f'{x}, del {find_nearest(nearest_num,x)} from Pump2')
  Pump_all_no_init = [*Pump1,*Pump2]
  Pump_all_no_init = np.array(sorted(Pump_all_no_init))
  return Pump1,Pump2,Pump_all_no_init



def find_match_isx_GPIO_csv(file_path,GPIO_folder=None):
    mouse_id = file_path.split('_Inscopix-')[0].split('/')[-1]
    date     = file_path.split('_Inscopix-')[1].split('_')[0]
    date_nodot = date.replace(".", "")
    if GPIO_folder == None:
      GPIO_folder = '/Users/hunglo/Documents/inscopix_csv/Ca2_csv/'
    # print(GPIO_folder)
    subfolders = sorted([x[0] for x in os.walk(GPIO_folder)])
    subfolders = [f for f in subfolders if mouse_id in f][0]
    # print(subfolders)
    extension = 'csv'
    os.chdir(subfolders) # cd to mouse folder
    result = sorted(glob.glob('*.{}'.format(extension)))
    result = [r for r in result if date_nodot in r and 'GPIO' in r][0]
    file_path_GPIO = os.path.join(subfolders,result)
    return file_path_GPIO

def findandload_match_isx_csv(file_path, time_difference, isx_folder = None):
    mouse_id = file_path.split('_Inscopix-')[0].split('/')[-1]
    date     = file_path.split('_Inscopix-')[1].split('_')[0]
    date_nodot = date.replace(".", "")
    if isx_folder == None:
      isx_folder = '/Users/hunglo/Documents/inscopix_csv/Ca2_csv/'
    subfolders = sorted([x[0] for x in os.walk(isx_folder)])
    subfolders = [f for f in subfolders if mouse_id in f][0]
    # print(subfolders)
    extension = 'csv'
    os.chdir(subfolders) # cd to mouse folder
    result = sorted(glob.glob('*.{}'.format(extension)))
    result = [r for r in result if date_nodot in r and '_celltraces.' in r]

    if 'BES0224' in mouse_id:
        # print('here check')
        df_merged = pd.DataFrame([])
        for patch in result:
            patch_num = '_'+str(patch.split('_')[1])
            df_all = pd.read_csv(subfolders+'/'+patch,header=[0,1],index_col=0)
            # Selecet only the accepted cells
            df_accepted = df_all.xs(' accepted',level='Time(s)/Cell Status',axis=1)
            df_accepted.columns = [str(col) + patch_num for col in df_accepted.columns]
            df_merged[df_accepted.columns] = df_accepted
        df_accepted = pd.DataFrame([])
        df_accepted = df_merged.copy()
    elif  'DSC016006' in mouse_id:
        # print('here check2')
        df_merged = pd.DataFrame([])
        for patch in result:
            patch_num = '_'+str(patch.split('_')[1])
            df_all = pd.read_csv(subfolders+'/'+patch,header=[0,1],index_col=0)
            # Selecet only the accepted cells
            df_accepted = df_all.xs(' accepted',level='Time(s)/Cell Status',axis=1)
            df_accepted.columns = [str(col) + patch_num for col in df_accepted.columns]
            df_merged[df_accepted.columns] = df_accepted
        df_accepted = pd.DataFrame([])
        df_accepted = df_merged.copy()

    else:
        df_all = pd.read_csv(os.path.join(subfolders,result[0]),header=[0,1],index_col=0)
        # Selecet only the accepted cells
        df_accepted = df_all.xs(' accepted',level='Time(s)/Cell Status',axis=1)
        if 'SNA' in mouse_id:
          df_accepted = df_accepted[15:] # drop first 15 seconds
          df_accepted = correct_baseline_modploy(df_accepted) # correct for baseline drift

    df_accepted.index = df_accepted.index+time_difference
    if time_difference>0:
      print('!!!!!!!!!!!!!!!!!!!!!!!!!!')
      print(f'time diff is positive (should be negative)!!!!!!!!!!!!!!!!!!!!!!\ncheck again manually:\n{mouse_id}-{date}')
      print('!!!!!!!!!!!!!!!!!!!!!!!!!!')
      if mouse_id in ['DSC012972','SNA089408']:
         pass # one recording actually is positive, but the values are all correct so skip
      else:
        sys.exit()

    return df_accepted


def load_GPIO_csv(file_path):
    file_path_GPIO = find_match_isx_GPIO_csv(file_path)
    
    if os.path.isfile(file_path_GPIO):
      ## here, check if the corresponding GPIO file exist
      df_GPIO = pd.read_csv(file_path_GPIO, header=[0], index_col=0)
    else:
       print(f"{file_path_GPIO.split('/')[-1]} doesn't exist, check again!")
    return df_GPIO

def extract_GPIO_trace(df_GPIO,plot=False):
    """
    Extract GPIO traces and binarized them to timestamps
    """
    GPIO_channel_list = [name for name in df_GPIO[' Channel Name'].unique() if 'IO' in name]
    ## this should extract both GPIO or IO prefix in difference nVista system
    gpio_pump_time_all = {}
    for channel in GPIO_channel_list:
        gpio_value = df_GPIO.loc[df_GPIO[' Channel Name'] == channel][' Value'].astype(float)
        # print(channel)
        ## check if there is signal at all
        if gpio_value.max()/(gpio_value.min()+1e-8) > 2: # prevent min is 0
            gpio_norm = (gpio_value-gpio_value.min())/(gpio_value.max()-gpio_value.min())
            gpio_timestamps = np.concatenate([np.zeros(1),gpio_norm[gpio_norm<0.5].index],axis=0)
            gpio_time_diff = np.diff(gpio_timestamps)
            gpio_pump_time = gpio_timestamps[1:][gpio_time_diff>0.04]
            gpio_pump_time_all[channel] = gpio_pump_time
            # print(f'First pump time: {gpio_pump_time_all[channel][0]}')

            if plot:
                plt.figure(figsize=[10,2])
                plt.title(channel)
                plt.eventplot(gpio_pump_time_all[channel],colors='C1',lineoffsets=-0.5,linewidths=0.5,label='pump activate timestamp')
                plt.plot(gpio_norm,lw=0.5,label='raw normalized')
                plt.legend(loc=4)
                plt.yticks([])
                plt.ylim([-1,1])
                sns.despine()
                plt.tight_layout()
                plt.xlim([320,340])
                # plt.xlim([410,420])
                plt.show()
        else:
            # print(f'No inputs in{channel} channel, skipped\n')
            # print(channel)
            pass
    ## Check if GPIO signals are missing, if so add an empty np array to the corresponding channel:
    # print(f'chan #{len(gpio_pump_time_all.keys())}')
    if len(gpio_pump_time_all.keys())==1:
        if ' GPIO-1' in gpio_pump_time_all.keys():
            gpio_pump_time_all[' GPIO-2'] = np.array([])
        elif ' IO1' in gpio_pump_time_all.keys():
            gpio_pump_time_all[' IO2'] = np.array([])
        elif ' IO2' in gpio_pump_time_all.keys():
            gpio_pump_time_all[' IO1'] = np.array([])
        elif ' GPIO-3' in gpio_pump_time_all.keys():
            gpio_pump_time_all[' GPIO-4'] = np.array([])
    ## simplify GPIO channel to just GPIO1 & 2
    keys = gpio_pump_time_all.keys()
    # print(f'chan #{len(gpio_pump_time_all.keys())}')
    # GPIO_pool = []
    for idx,key in enumerate(keys):
      # print(idx,key)
      # print(gpio_pump_time_all.get(key))
    #   GPIO_pool.append(np.min(gpio_pump_time_all[key]))
      if idx ==0:
        GPIO_pump1 = gpio_pump_time_all[key]
      else:
        GPIO_pump2 = gpio_pump_time_all[key]
    GPIO_pump1_clean, GPIO_pump2_clean, GPIO_pump_all_no_init = remove_isx_GPIO_init(GPIO_pump1,GPIO_pump2)
    # GPIO_min = np.min(GPIO_pump_all_no_init)
    # print(len(GPIO_pump1_clean), len(GPIO_pump2_clean))

    return GPIO_pump1_clean, GPIO_pump2_clean, GPIO_pump_all_no_init

def compare_pheno_inscopix_gpio(PAll,GPIO_pump_all_no_init):
    if len(PAll) <2:
        print('no enough drinking events')
        time_difference = 0
    else:
        pheno_min = np.min(PAll)
        GPIO_min  = np.min(GPIO_pump_all_no_init)
        time_difference = pheno_min-GPIO_min
        # print(pheno_min,GPIO_min)
        # print(f'pheno:{pheno_min}\nisx_GPIO:{GPIO_min}\ntime difference: {time_difference:.5f}')
    return time_difference

def get_time_difference(PAll,GPIO_pump_all_no_init,windowsize=100):
    """
    This is the current one I used now, use the pynapple cross-correlogram funciton, I can get the time difference without fully matching the timestamps
    """
    if len(PAll) <2:
      print('no enough drinking events')
      time_difference = 0
    else:
      ts1 = nap.Ts(np.array(PAll))
      ts2 = nap.Ts(np.array(GPIO_pump_all_no_init))
      ts1_time_array = ts1.index.values
      ts2_time_array = ts2.index.values

      # if date =='21.09.10':
      #    windowsize = 5

      binsize=0.001
      cc12, xt = nap.cross_correlogram(t1=ts1_time_array,t2=ts2_time_array,binsize=binsize,windowsize=windowsize)

      idx_max = np.argmax(cc12)
      print(f'time difference: {-xt[idx_max]}')
      time_difference = -xt[idx_max] # set to negative
      # plt.figure()
      # plt.bar(xt, cc12, binsize)
      # plt.xlabel("Time t1 (s)")
      # plt.ylabel("CC")
    return time_difference 

def loadandsync_incopix_csv(file_path,mouse_id):
    """
    Note here we use the phenosys csv file as the file_path, 
    and we should be able to directly load and sync corresponding Ca2+ csv files
    """

    # file_path_full = os.path.join(Phenosys_folder,file_path)
    _,_,_,PAll = load_phenosys(file_path)    
    df_GPIO = load_GPIO_csv(file_path)
    _,_, GPIO_pump_all_no_init = extract_GPIO_trace(df_GPIO)
    # time_difference = compare_pheno_inscopix_gpio(PAll,GPIO_pump_all_no_init)
    time_difference = get_time_difference(PAll,GPIO_pump_all_no_init)
    if abs(round(time_difference,0)) == 100:
      ## getting boundary cases, might due to too few events or too large window
      ## re-perform get time difference with small time window
      time_difference = get_time_difference(PAll,GPIO_pump_all_no_init,windowsize=20)
    elif '21.09.10' in file_path: # SNA089408 has some issues with x-corr, multiple max corr values, used a old function here
       print('unstable x-corr, used old time difference function')
       time_difference = compare_pheno_inscopix_gpio(PAll,GPIO_pump_all_no_init) # used the old one 
       print(f'time difference (corrected): {time_difference}')

    df_accepted = findandload_match_isx_csv(file_path, time_difference)
    print(f'==================\nnum of neurons: {df_accepted.shape[1]}')

    from scipy import stats
    df_z = df_accepted.apply(stats.zscore,axis=0)
    return df_accepted, df_z, time_difference

def correct_baseline_modploy(df_accepted):
    input_array = df_accepted.values
    polynomial_degree=2 #only needed for Modpoly and IModPoly algorithm
    baseObj=BaselineRemoval(input_array)
    Modpoly_output=baseObj.ModPoly(polynomial_degree)
    df_corrected = pd.DataFrame(data=Modpoly_output,index=df_accepted.index,columns=[' C1'])
    return df_corrected

def altspace(start, step, count, endpoint=False, **kwargs):
   stop = start+(step*count)
   return np.linspace(start, stop, count, endpoint=endpoint, **kwargs)

#@title Plotting setting 
## Plotting setting
# plt.rcParams["font.family"] = "Arial" # The default 'Dejavu' is fine
plt.rcParams.update({'font.size': 10})
my_color_map = ['#56b4e9',
                '#e69f00',
                '#009e73',
                '#f0e442',
                '#0072b2',
                '#d55e00',
                '#cc79a7']
# if show_color_map:
#     plt.figure(figsize=[5,1.5])
#     for idx, c in enumerate(my_color_map):
#         plt.plot(idx,0,color=c,marker='.',markersize=40)
#         plt.text(x=idx,y=1,s=idx,color='k',ha='center',va='center',size=20)
#     plt.title('custom color map')
#     plt.axis('off')
#     plt.ylim([-0.5,1.5])
#     plt.tight_layout()
#     plt.show()   

## For getting hex code from rgb
# rgb_code = [204,121,167]
# hex_color = '#%02x%02x%02x' % (rgb_code[0], rgb_code[1], rgb_code[2])
# print(hex_color)

#@title Functions for plotting

def title_plot_name(mouse_id, Food_deprivation, date):

  if Food_deprivation == 'severe':
      title = mouse_id+' '+date+' '+Genotype+ ' Neurons (severe FD)'
  elif Food_deprivation == 'mild':
      title = mouse_id+' '+date+' '+Genotype+ ' Neurons (mild FD)'
  elif Food_deprivation == 'none':
      title = mouse_id+' '+date+' '+Genotype+ ' Neurons ($\it{ad}$ $\it{libitum}$)'
  plot_name = mouse_id+'_'+date+'_'+Genotype+'_'
  return title, plot_name

def length_delivery(Pump1,Pump2):
  deliver_milk = len(Pump1)
  deliver_water = len(Pump2)
  return deliver_milk, deliver_water

def GPIO_eventplot(Pump1, Pump2, ax):
  if len(Pump2) == 0:
    GPIO_plot = [Pump1]
    ax.eventplot(GPIO_plot,color=['C3'],lineoffsets=[y_offset*1.1],linewidths=0.5, linelength=1,linestyle = 'None',alpha=1)
  else:
    GPIO_plot = [Pump1,Pump2]
    ax.eventplot(GPIO_plot,color=[my_color_map[2],my_color_map[1]],lineoffsets=[1,0],linewidths=0.5, linelength=0.8,linestyle = 'None',alpha=1)
  return ax

def Ca2_eventplot(df_event_sorted,ax):
  for i,name in enumerate(df_event_sorted.columns):
    data=df_event_sorted[name][df_event_sorted[name] > 0]
    ax.eventplot(data.index.values,lineoffsets = i+1,linewidth=0.3,linelength=1,alpha=1,color='C0')#my_color_map[0])
  return ax

####



############# testing code

# Phenosys_folder = '/Users/hunglo/Documents/inscopix_csv/Phenosys_csv/'
# extension = 'csv'
# os.chdir(Phenosys_folder)
# result = sorted(glob.glob('*.{}'.format(extension)))
# # result = [r for r in result if '21.01.06' in r]
# # result = [r for r in result if 'DSC012849' in r]

# for file_path in tqdm(result):
#   mouse_id = file_path.split('_Inscopix-')[0].split('/')[-1]
#   date     = file_path.split('_Inscopix-')[1].split('_')[0]
#   print(mouse_id,date)
#   L1,P1,P2,PAll = load_phenosys(file_path)

#   L1_,P1_,P2_,PAll_ = load_phenosys_nofilter(file_path)

#   if 0:
#     x_shift = PAll_[-1]
#     if x_shift< 1550:
#        x_shift = 1550 # force output to be 8 panels for consistency
#     x_bar = [x*180+300 for x in range(18)]
#     x_bar = [0] + x_bar
#     x_bar = np.array(x_bar)
#     x_bar = x_bar[x_bar<x_shift]
#     n=0
#     plt.figure(figsize=[12,4])
#     plt.suptitle(f'{mouse_id}_{date}')
#     for x in x_bar:
#       ax = plt.subplot(2, (len(x_bar)+1)//2, n+1)
#       ax.eventplot([P1,P2,P1_,P2_],linewidths=1,colors=['C0','C0','C1','C1'],linelengths=0.8)
#       ax.set_xlim([x-10,x+130])
#       sns.despine(ax=ax)
#       ax.set_yticks([])
#       n+=1
#     plt.tight_layout()
#     plt.show()
#   # plt.figure()
#   # plt.eventplot([P2,P2_],linewidths=0.8,colors=['C1','C1'])
#   # # plt.show()

#   # print(len(P1),len(P1_))
#   # print(len(P2),len(P2_))

